###############################################################################
# NZ POLL FORECAST MODEL
#
# Workflow:
#   1. Load and prepare polling data
#   2. Fit weighted LOESS trends
#   3. Estimate historical polling covariance
#   4. Run Monte Carlo election simulations
#   5. Summarise coalition probabilities
#
###############################################################################
###############################################################################
# MODEL ASSUMPTIONS
#
# 1. Polling residuals are representative of future polling uncertainty.
# 2. Residuals follow an approximately multivariate normal distribution.
# 3. LOESS trends are a reasonable estimate of current voting intention.
# 4. Poll sample size is an appropriate weighting variable.
# 5. Māori electorate outcomes can be approximated from TPM party vote.
# 6. Coalition classifications reflect predefined political assumptions.
# 7. No explicit house effects, turnout effects, or late campaign swings are
#    modelled.
# 8. Party Electorate Distribution, which can can effect overhang, is 
#    reasonably accurate (manually set line 896)
#
###############################################################################

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(MASS)

###############################################################################
# DATA PREPARATION
#
# Purpose:
#   Read polling observations and standardise them for modelling.
#
# Notes:
#   - Midate represents the midpoint fieldwork date for each poll and is used
#     as the primary time variable throughout the model.
#   - A numeric day counter is created because LOESS requires a continuous
#     predictor rather than a Date object.
#   - Polls are sorted chronologically to ensure trend estimation and the
#     selection of the most recent trend value occur in time order.
#   - The first row is being removed as the previous election included in the
#     wikipedia polling for the general election page I got the data from
#
# Output:
#   A polling dataset ordered by date with a numeric time index ('days')
#   suitable for trend fitting.
###############################################################################

prepare_polls <- function(
    path,
    remove_first_row = TRUE
) {
  
  polls <- read_csv(path)
  
  if (remove_first_row) {
    polls <- polls[-1, ]
  }
  
  polls |>
    mutate(
      Midate = dmy(Midate)
    ) |>
    arrange(Midate) |>
    mutate(
      `Polling organisation` = trimws(`Polling organisation`),
      days = as.numeric(
        Midate - min(Midate)
      )
    )
  
}

###############################################################################
# LOESS TRENDS
#
# Purpose:
#   Estimate the underlying level of support for each party while smoothing
#   out poll-to-poll noise.
#
# Method:
#   - A separate LOESS model is fitted for each party.
#   - Poll sample size is used as a weight so larger polls exert greater
#     influence on the fitted trend than smaller polls.
#   - The span parameter controls the amount of smoothing:
#       * Smaller values follow short-term movement more closely.
#       * Larger values produce smoother, less responsive trends.
#
# Why this matters:
#   Raw polling contains sampling variation and methodological differences.
#   The LOESS trend is treated as the best estimate of current support and
#   becomes the centre point of the election simulations.
#
# Outputs:
#   - Polling dataset with a trend column for each party.
#   - latest_poll vector containing the most recent trend estimate for every
#     party, which becomes the simulation baseline.
###############################################################################

fit_poll_trends <- function(
    polls,
    parties,
    span = 0.30,
    weight_col = "Sample size"
) {
  
  latest_poll <- numeric(length(parties))
  names(latest_poll) <- parties
  
  for (p in parties) {
    
    fit <- loess(
      stats::as.formula(
        paste0(p, " ~ days")
      ),
      data = polls,
      weights = polls[[weight_col]],
      span = span
    )
    
    trend_col <- paste0(p, "_trend")
    
    polls[[trend_col]] <- predict(fit)
    
    latest_poll[p] <- tail(
      polls[[trend_col]],
      1
    )
  }
  
  list(
    polls = polls,
    latest_poll = latest_poll
  )
  
}

###############################################################################
# RAW POLL RESIDUALS
#
# Purpose:
# Measure deviation from the LOESS trend before adjusting for pollster bias.
#
# Definition:
#
# Raw Residual =
# Observed Poll
# - LOESS Trend
#
# Components contained within residuals:
#
# - Sampling error
# - Campaign noise
# - Pollster house effects
# - Unexplained variation
###############################################################################

add_poll_residuals <- function(
    polls,
    parties
) {
  
  for (p in parties) {
    
    polls[[paste0(p, "_resid")]] <-
      polls[[p]] -
      polls[[paste0(p, "_trend")]]
  }
  
  polls
  
}

###############################################################################
# HOUSE EFFECT ESTIMATION
#
# Purpose:
#   Estimate systematic pollster-specific bias relative to trend.
#
# Method:
#   - Calculate average residual for each polling organisation.
#   - Apply empirical-Bayes style shrinkage toward zero.
#   - Larger poll histories receive less shrinkage.
#
# Shrinkage:
#
#   adjusted =
#       raw_effect *
#       n / (n + shrinkage_k)
#
# Output:
#   Party-specific house effect columns.
###############################################################################

estimate_house_effects <- function(
    polls,
    parties,
    organisation_col = "Polling organisation",
    weight_col = "Sample size",
    shrinkage_k = 10
) {
  
  house_effect_tables <- list()
  
  for (p in parties) {
    
    resid_col <- paste0(p, "_resid")
    
    house_table <-
      polls |>
      group_by(.data[[organisation_col]]) |>
      summarise(
        n = n(),
        
        raw_house =
          weighted.mean(
            .data[[resid_col]],
            .data[[weight_col]],
            na.rm = TRUE
          ),
        
        .groups = "drop"
      ) |>
      
      mutate(
        house_effect =
          raw_house *
          n / (n + shrinkage_k)
      )
    
    house_effect_tables[[p]] <- house_table
    
    effect_lookup <-
      setNames(
        house_table$house_effect,
        house_table[[organisation_col]]
      )
    
    polls[[paste0(p, "_house")]] <-
      effect_lookup[
        polls[[organisation_col]]
      ]
  }
  
  list(
    polls = polls,
    house_effect_tables = house_effect_tables
  )
}

###############################################################################
# HOUSE-ADJUSTED RESIDUALS
#
# Purpose:
#   Remove systematic pollster effects from residuals.
#
# Components remaining:
#
#   - Sampling error
#   - Campaign noise
#   - Unexplained variation
#
# These adjusted residuals form the basis of covariance estimation.
###############################################################################

add_adjusted_residuals <- function(
    polls,
    parties
) {
  
  for (p in parties) {
    
    polls[[paste0(p, "_adj_resid")]] <-
      polls[[paste0(p, "_resid")]] -
      polls[[paste0(p, "_house")]]
  }
  
  polls
}


###############################################################################
# COVARIANCE ESTIMATION
#
# Purpose:
#   Estimate the historical joint uncertainty structure of polling errors.
#
# Why covariance matters:
#   Election outcomes depend not only on uncertainty for individual parties,
#   but also on how party support moves together.
#
# Examples:
#   - If Labour overperforms trend, Green may also overperform trend.
#   - If National gains support, another centre-right party may lose support.
#
# Method:
#   - Only polls with complete residual information are retained.
#   - Poll sample size is used as a weight so larger polls contribute more to
#     the covariance estimate.
#   - cov.wt() produces a maximum-likelihood covariance matrix describing
#     historical residual behaviour.
#
# Outputs:
#   covariance:
#       Joint variance-covariance matrix used directly in simulations.
#
#   correlation:
#       Standardised version used for diagnostics and interpretation.
#
#   residual_data:
#       Retained to support validation and diagnostic checks.
###############################################################################


estimate_poll_covariance <- function(
    polls,
    parties,
    weight_col = "Sample size"
) {
  
  residual_cols <-
    paste0(
      parties,
      "_adj_resid"
    )
  
  complete_cases <-
    complete.cases(
      polls[, residual_cols]
    )
  
  residual_data <-
    polls[
      complete_cases,
      residual_cols
    ]
  
  weights <-
    polls[[weight_col]][
      complete_cases
    ]
  
  covariance <-
    cov.wt(
      x = residual_data,
      wt = weights,
      method = "ML"
    )$cov
  
  list(
    covariance = covariance,
    correlation = cov2cor(
      covariance
    ),
    residual_data = residual_data
  )
}

###############################################################################
# COVARIANCE DIAGNOSTICS
#
# Purpose:
#   Verify that the covariance matrix is suitable for multivariate simulation.
#
# Diagnostics:
#   Eigenvalues:
#       Used to assess whether the covariance matrix is positive definite.
#
#   Positive-definite test:
#       Required by multivariate normal simulation.
#       Negative or zero eigenvalues indicate mathematical problems that may
#       cause simulation failures or unrealistic behaviour.
#
#   Residual standard deviations:
#       Show the magnitude of historical polling error for each party.
#
# Why this matters:
#   A covariance matrix that fails these checks can invalidate the Monte Carlo
#   model and lead to unstable simulation results.
###############################################################################

covariance_diagnostics <- function(
    covariance,
    residual_data
) {
  
  eigenvalues <-
    eigen(covariance)$values
  
  residual_sd <-
    sapply(
      residual_data,
      sd,
      na.rm = TRUE
    )
  
  list(
    eigenvalues = eigenvalues,
    positive_definite = all(
      eigenvalues > 0
    ),
    residual_sd = residual_sd
  )
  
}

###############################################################################
# HOUSE EFFECT DIAGNOSTICS
#
# Purpose:
#   Display estimated pollster effects by party.
###############################################################################

display_house_effects <- function(
    house_effect_tables,
    parties
) {
  
  for (p in parties) {
    
    cat("\n")
    cat(
      "================================================\n"
    )
    
    cat(
      paste(
        "HOUSE EFFECTS:",
        p
      )
    )
    
    cat("\n")
    cat(
      "================================================\n"
    )
    
    print(
      house_effect_tables[[p]] |>
        arrange(
          desc(abs(house_effect))
        )
    )
  }
}

###############################################################################
# TREND PLOTS
#
# Purpose:
#   Provide a visual validation of the LOESS trend model.
#
# Interpretation:
#   Grey line:
#       Observed polling results.
#
#   Blue line:
#       Smoothed estimate of underlying support.
#
# Analysts should review these plots to confirm:
#   - Trends are not over-smoothed or under-smoothed.
#   - Major political events are reasonably reflected.
#   - End-of-series behaviour appears credible.
#
# Output:
#   Named list of ggplot objects, one for each party.
###############################################################################

create_trend_plots <- function(
    polls,
    parties
) {
  
  lapply(
    parties,
    function(p) {
      
      ggplot(
        polls,
        aes(x = Midate)
      ) +
        
        geom_line(
          aes(
            y = .data[[p]]
          ),
          colour = "grey70"
        ) +
        
        geom_line(
          aes(
            y = .data[[paste0(
              p,
              "_trend"
            )]]
          ),
          colour = "blue",
          linewidth = 1
        ) +
        
        theme_minimal() +
        
        labs(
          title = paste(
            p,
            "LOESS Trend"
          ),
          y = "Percent"
        )
      
    }
  ) |>
    setNames(parties)
  
}

###############################################################################
# VOTE SIMULATION
#
# Purpose:
#   Generate one plausible election-day vote distribution.
#
# Method:
#   - The latest LOESS trend forms the expected vote share.
#   - Historical polling covariance defines uncertainty around that estimate.
#   - MASS::mvrnorm() draws a random outcome from the joint distribution.
#
# Assumptions:
#   - Polling error follows an approximately multivariate normal distribution.
#   - Historical polling behaviour remains informative of future uncertainty.
#
# Post-processing:
#   - Negative simulated vote shares are truncated to zero.
#   - Vote shares are renormalised to sum to 100%.
#
# Output:
#   A complete simulated party vote vector.
###############################################################################

simulate_vote <- function(
    latest_poll,
    covariance
) {
  
  vote <-
    MASS::mvrnorm(
      n = 1,
      mu = latest_poll,
      Sigma = covariance
    )
  
  names(vote) <- names(latest_poll)
  
  vote[vote < 0] <- 0
  
  vote / sum(vote) * 100
  
}

###############################################################################
# MAORI ELECTORATE MODEL
#
# This is the bit I have least confidence in representing the actual world
# as it is both a small number of electorates, different to the others, and
# poorly canvassed. So this is more of a guesstimate.
#
# Purpose:
#   Model uncertainty in Māori electorate outcomes.
#
# Assumption:
#   Te Pāti Māori electorate performance is linked to party vote strength.
#
# Mechanism:
#   - TPM vote share is converted into a win probability using a logistic
#     function.
#   - Remaining Māori seats are then simulated using a binomial draw.
#
# Notes:
#   - This is a simplified behavioural model rather than a seat-by-seat
#     electorate forecast.
#   - The coefficients within the logistic equation are modelling assumptions
#     and should ideally be calibrated using historical electorate results.
#
# Output:
#   Seat counts attributed to Labour, TPM, and Independent candidates.
###############################################################################

simulate_maori_electorates <- function(
    tpm_vote,
    independent_count
) {
  
  remaining <- 7 - independent_count
  
  probability <-
    plogis(
      -0.5 + 0.4 * tpm_vote
    )
  
  tpm_seats <-
    rbinom(
      1,
      remaining,
      probability
    )
  
  c(
    LAB = remaining - tpm_seats,
    TPM = tpm_seats,
    IND = independent_count
  )
  
}

###############################################################################
# SAINTE-LAGUE ALLOCATION
#
# Purpose:
#   Convert party vote estimates into parliamentary seat allocations under
#   New Zealand's MMP system.
#
# Qualification rules:
#   Parties qualify if they:
#       - Reach the 5% party vote threshold; or
#       - Win at least one electorate seat.
#
# Method:
#   - Generate Sainte-Laguë quotients for all qualifying parties.
#   - Rank all quotients nationally.
#   - Allocate seats to the highest quotients until Parliament is filled.
#
# Overhang handling:
#   If a party wins more electorate seats than its proportional entitlement,
#   those electorate seats are retained and Parliament expands accordingly.
#
# Outputs:
#   - Final seat allocation.
#   - Number of overhang seats.
#   - Effective Parliament size.
###############################################################################

allocate_seats <- function(
    vote_share,
    electorate_seats,
    parliament_size = 120
) {
  
  parties <- names(electorate_seats)
  
  qualifying <-
    vote_share >= 5 |
    electorate_seats > 0
  
  quotients <- list()
  
  for (p in names(vote_share)) {
    
    if (!qualifying[p]) next
    
    divisors <- c(
      1.4,
      seq(3, 299, by = 2)
    )
    
    quotients[[p]] <-
      data.frame(
        party = p,
        quotient =
          vote_share[p] / divisors
      )
  }
  
  qtab <- bind_rows(quotients)
  
  winners <-
    qtab |>
    arrange(desc(quotient)) |>
    slice_head(
      n = parliament_size
    )
  
  proportional <-
    table(winners$party)
  
  seat_alloc <-
    setNames(
      rep(0, length(parties)),
      parties
    )
  
  seat_alloc[
    names(proportional)
  ] <- as.vector(proportional)
  
  seats <- seat_alloc
  overhang <- 0
  
  for (p in parties) {
    
    if (electorate_seats[p] > seats[p]) {
      
      extra <-
        electorate_seats[p] -
        seats[p]
      
      overhang <-
        overhang + extra
      
      seats[p] <-
        electorate_seats[p]
    }
    
  }
  
  list(
    seats = seats,
    overhang = overhang,
    parliament = parliament_size + overhang
  )
  
}

###############################################################################
# GOVERNMENT CLASSIFICATION
#
# Purpose:
#   Translate seat allocations into politically meaningful government outcomes.
#
# Coalition assumptions:
#   Right bloc:
#       National + ACT + NZ First
#
#   Left bloc:
#       Labour + Green + TPM + Independents
#       (this assumes a minority government w. support may be possible)
#
#   TOP:
#       Treated as a potential cross-bloc coalition partner.
#
# Logic:
#   - Check whether either bloc already holds a majority.
#   - If not, determine whether TOP can create a majority.
#   - If both blocs could govern with TOP, classify as "TOP Decides".
#   - Otherwise classify as "Hung".
#
###############################################################################

classify_government <- function(
    seats,
    parliament
) {
  
  majority <-
    floor(parliament / 2) + 1
  
  right_bloc <-
    sum(
      seats[c(
        "NAT",
        "ACT",
        "NZF"
      )]
    )
  
  left_bloc <-
    sum(
      seats[c(
        "LAB",
        "GRN",
        "TPM",
        "IND"
      )]
    )
  
  top_seats <-
    as.numeric(
      seats["TOP"]
    )
  
  if (right_bloc >= majority) {
    return("National Majority")
  }
  
  if (left_bloc >= majority) {
    return("Labour Majority")
  }
  
  if (top_seats == 0) {
    return("Hung")
  }
  
  right_possible <-
    right_bloc + top_seats >= majority
  
  left_possible <-
    left_bloc + top_seats >= majority
  
  if (right_possible & left_possible) {
    return("TOP Decides")
  }
  
  if (right_possible) {
    return("National + TOP")
  }
  
  if (left_possible) {
    return("Labour + TOP")
  }
  
  "Hung"
  
}

###############################################################################
# SINGLE SIMULATION
#
# Purpose:
#   Execute one end-to-end election scenario.
#
# Workflow:
#   1. Simulate a national (small n) party vote.
#   2. Simulate Māori electorate outcomes.
#   3. Update electorate seat counts.
#   4. Allocate parliamentary seats.
#   5. Classify the resulting government formation.
#
# Why this function exists:
#   It combines all individual model components into a single contestable
#   election outcome which can be repeated thousands of times.
#
# Output:
#   A complete simulated parliament and government classification.
###############################################################################

simulate_once <- function(
    latest_poll,
    covariance,
    base_electorates,
    independent_count
) {
  
  vote <-
    simulate_vote(
      latest_poll,
      covariance
    )
  
  electorates <-
    base_electorates
  
  maori <-
    simulate_maori_electorates(
      vote["TPM"],
      independent_count
    )
  
  electorates["LAB"] <-
    electorates["LAB"] +
    maori["LAB"]
  
  electorates["TPM"] <-
    maori["TPM"]
  
  electorates["IND"] <-
    maori["IND"]
  
  allocation <-
    allocate_seats(
      vote_share = c(
        vote,
        IND = 0
      ),
      electorate_seats =
        electorates
    )
  
  list(
    seats = allocation$seats,
    parliament = allocation$parliament,
    overhang = allocation$overhang,
    outcome = classify_government(
      allocation$seats,
      allocation$parliament
    )
  )
  
}

###############################################################################
# SCENARIO SIMULATION
#
# Purpose:
#   Estimate probabilities of different government outcomes.
#
# Method:
#   - Repeat the full election simulation many times.
#   - Record the resulting government classification.
#   - Convert counts into probabilities.
#
# Interpretation:
#   If "National Majority" appears in 4,500 of 10,000 simulations,
#   the estimated probability is 45%.
#
# n_sims:
#   Larger values improve stability of probability estimates but increase
#   processing time.
#
# Output:
#   Probability distribution across government formation scenarios.
###############################################################################

run_scenario <- function(
    latest_poll,
    covariance,
    base_electorates,
    independent_count,
    n_sims = 10000
) {
  
  simulations <-
    replicate(
      n_sims,
      simulate_once(
        latest_poll = latest_poll,
        covariance = covariance,
        base_electorates = base_electorates,
        independent_count = independent_count
      ),
      simplify = FALSE
    )
  
  outcomes <-
    table(
      vapply(
        simulations,
        `[[`,
        character(1),
        "outcome"
      )
    )
  
  outcomes / sum(outcomes)
  
}

###############################################################################
# SCENARIO REPORTING
#
# Purpose:
#   Present simulation probabilities in a readable console format.
#
# Design:
#   - Ensures all standard outcome categories are displayed.
#   - Missing categories are reported as 0%.
#   - Consistent formatting supports comparison across scenarios.
#
# Typical use:
#   Compare how government probabilities change under differing assumptions
#   about the number of independent Māori electorate MPs.
###############################################################################

display_scenario <- function(
    results,
    independent_count
) {
  
  categories <- c(
    "National Majority",
    "National + TOP",
    "Labour Majority",
    "Labour + TOP",
    "TOP Decides",
    "Hung"
  )
  
  cat("\n")
  cat("====================================================\n")
  cat(
    paste(
      "SCENARIO:",
      independent_count,
      "INDEPENDENT MAORI MP(S)"
    )
  )
  cat("\n")
  cat("====================================================\n")
  
  for (k in categories) {
    
    p <- ifelse(
      k %in% names(results),
      results[[k]],
      0
    )
    
    cat(
      sprintf(
        "%-20s %6.2f%%\n",
        k,
        100 * p
      )
    )
  }
  
}

###############################################################################
# ORCHESTRATION
#
# Purpose:
#   Execute the complete forecasting pipeline from polling data through to
#   election simulations.
#
# Pipeline overview:
#   1. Define modelled parties and electorate assumptions.
#   2. Import and prepare polling data.
#   3. Estimate party trends.
#   4. Calculate polling residuals.
#   5. Estimate covariance structure.
#   6. Validate covariance quality.
#   7. Produce diagnostics and visualisations.
#   8. Export covariance outputs.
#   9. Run election simulations under multiple Māori electorate scenarios.
#  10. Report government probabilities.
#
# Key modelling concept:
#   The latest smoothed polling trend provides the expected election outcome,
#   while historical polling residuals provide the uncertainty around that
#   expectation.
###############################################################################

parties <- c(
  "NAT",
  "LAB",
  "GRN",
  "ACT",
  "NZF",
  "TPM",
  "TOP"
)

# Lookout, this method is really sensitive to expected electorate results
# I just looked a the Lab vs Nat polling and picking some arbitrary numbers

base_electorates <- c(
  NAT = 35,
  LAB = 28,
  GRN = 3,
  ACT = 1,
  NZF = 0,
  TPM = 3,
  TOP = 0,
  IND = 1
)

polls <-
  prepare_polls(
    path = "allpolls.csv"
  )

trend_results <-
  fit_poll_trends(
    polls,
    parties,
    span = 0.30
  )

polls <-
  trend_results$polls

latest_poll <-
  trend_results$latest_poll

polls <-
  add_poll_residuals(
    polls,
    parties
  )

house_results <-
  estimate_house_effects(
    polls,
    parties,
    shrinkage_k = 10
  )

polls <-
  house_results$polls

polls <-
  add_adjusted_residuals(
    polls,
    parties
  )

cov_results <-
  estimate_poll_covariance(
    polls,
    parties
  )

cov_matrix <-
  cov_results$covariance

cor_matrix <-
  cov_results$correlation

diagnostics <-
  covariance_diagnostics(
    covariance = cov_matrix,
    residual_data =
      cov_results$residual_data
  )

display_house_effects(
  house_results$house_effect_tables,
  parties
)


###############################################################################
# DIAGNOSTIC OUTPUT
#
# Purpose:
# Display key model diagnostics for validation and review.
#
# Outputs include:
# - Current trend estimates.
# - Polling covariance and correlation matrices.
# - Covariance eigenvalues.
# - Historical residual standard deviations.
###############################################################################

print(round(latest_poll, 2))

print(round(cov_matrix, 4))

print(round(cor_matrix, 2))

print(round(
  diagnostics$eigenvalues,
  6
))

print(round(
  diagnostics$residual_sd,
  3
))

###############################################################################
# PLOTS
#
# Purpose:
# Generate and display trend charts for each modelled party.
#
# These plots provide a quick visual check of observed polling against the
# fitted LOESS trends.
###############################################################################

trend_plots <-
  create_trend_plots(
    polls,
    parties
  )

invisible(
  lapply(
    trend_plots,
    print
  )
)

###############################################################################
# OUTPUT FILES
#
# Purpose:
# Save covariance and correlation matrices for reference, validation, or
# reuse in later analysis.
###############################################################################

write.csv(
  cov_matrix,
  "covariance_matrix.csv"
)

write.csv(
  cor_matrix,
  "correlation_matrix.csv"
)

###############################################################################
# SIMULATIONS
#
# Purpose:
# Run election simulations under alternative assumptions regarding the
# number of independent Māori electorate MPs.
#
# Results are summarised as probabilities of each government outcome.
###############################################################################

scenario_0 <-
  run_scenario(
    latest_poll,        
    cov_matrix,
    base_electorates,
    independent_count = 0
  )

scenario_1 <-
  run_scenario(
    latest_poll,
    cov_matrix,
    base_electorates,
    independent_count = 1
  )

scenario_2 <-
  run_scenario(
    latest_poll,
    cov_matrix,
    base_electorates,
    independent_count = 2
  )

display_scenario(
  scenario_0,
  0
)

display_scenario(
  scenario_1,
  1
)

display_scenario(
  scenario_2,
  2
)

