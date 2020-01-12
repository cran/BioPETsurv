\name{sim_data}
\alias{sim_data}

\title{Simulating Biomarker and Survival Observations}
\description{
This function simulates biomarkers and generates survival observations depending on biomarker values. The simulated data can be used to explore prognostic enrichment using \code{surv_enrichment}.
}
\usage{
sim_data(n = 500, biomarker = "normal", effect.size = 1.25,
         baseline.hazard = "constant", end.time = 10,
         end.survival = 0.5, shape = NULL, seed = 2333)
}
\arguments{
  \item{n}{
  The number of observations to simulate.
}
  \item{biomarker}{
  Character specifying the shape of the biomarker distribution. Choices are \code{normal} for a symmetric distribution and \code{lognormal} for a right-skewed distribution.
}
  \item{effect.size}{
  The hazard ratio corresponding to one standard deviation increment in the biomarker.
}
  \item{baseline.hazard}{
  Character ("constant"/"increasing"/"decreasing") specifying whether the overall hazard in the population is constant, increasing or decreasing over time.
}
  \item{end.time}{
  The length of observation in the simulated dataset. In the data simulation, any events after this time will be censored at this time.
}
  \item{end.survival}{
  The survival rate in the population at the end of observation.
}
  \item{shape}{
  (Optional) the Weibull shape parameter for the baseline hazard. Values smaller and larger than 1 correspond to decreasing and increasing respectively.
}
  \item{seed}{
  (Optional) specify the random seed used for simulation.
}
}
\details{
  The biomarker will be simulated from a standardized normal or lognormal distribution. It is important that \code{effect.size} should correspond to a 1 SD increment in the biomarker. Conditioning on the biomarker values and assuming proportional hazards, survival times are simulated from a Weibull distribution with user-specified shape parameter, and the scale parameter is determined by the specified event rate and effect size.
}
\value{
  Returns a list of the following items:
  \item{data}{
    A data frame with 4 columns: the value of biomarker, observed event time, event indicator and the true event time.
  }
  \item{km.plot}{
    The Kaplan-Meier survival curves of the simulated dataset at enrichment levels 0, 25\%, 50\% and 75\%.
  }
}
\examples{
  ## Simulate a dataset with 500 observations,
  ## where the biomarker is Normally distributed (with SD=1).
  ## The hazard ratio corresponding to every one unit of increament in the biomarker is 1.25.
  ## The observation period is 10 months,
  ## and the survival probability of the population at the end of observation is 0.5.
  ## Hazards are constant over time.
  sim_obj <- sim_data(n = 500, biomarker = "normal", effect.size = 1.25,
                     baseline.hazard = "constant", end.time = 10, end.survival = 0.5)
  dat <- sim_obj$data
}
