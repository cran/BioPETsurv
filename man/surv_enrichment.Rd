\name{surv_enrichment}
\alias{surv_enrichment}
\title{Prognostic Enrichment for Clinical Trials with Survival Outcomes}
\description{
This function evaluates biomarkers for prognostic enrichment of clinical trials with survival outcomes, using data with biomarker values and survival observations.
A more detailed tutorial is available \href{https://github.com/chengs94/BioPETsurv/blob/master/README.md}{here}.
}
\usage{
surv_enrichment(formula, data, hr = 0.8, end.of.trial=NULL, a=NULL, f=NULL,
     method = "KM", lambda = 0.05,
     cost.screening = NULL, cost.keeping = NULL, cost.unit.keeping = NULL,
     power = 0.9, alpha = 0.05, one.sided = FALSE,
     selected.biomarker.quantiles = seq(from = 0, to = 0.95, by = 0.05),
     do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333,
     print.summary.tables = FALSE)
}
\arguments{
  \item{formula}{
    Object of class \code{formula}, in the form \code{outcome} ~ \code{predictors}, where the outcome is a survival object as returned by function \code{Surv()}. The predictor(s) include the biomarker(s) of interest and/or other covariates. If multiple predictors are inputed, this function will consider a composite biomarker based on a Cox model with all predictors.
  }
  \item{data}{
    Data frame containing the survival outcome (as returned by \code{Surv()}) and predictors specified in the \code{formula} argument. Observations with missing value(s) will be dropped.
  }
  \item{hr}{
    The hazard ratio (comparing the treatment versus control group) that the trial seeks to detect. Should be a number between 0 and 1.
  }
  \item{end.of.trial}{
    A scalar or a vector that gives the duration(s) of fixed-length trial(s) that are being considered. The unit of time should be the same as the survival times in the data. Either \code{end.of.trial} or a combination of \code{a} and \code{f} must be specified.
  }
  \item{a}{
    A scalar specifying the accrual period of the trial. It is assumed that patients are recruited at a uniform rate during this period.
  }
  \item{f}{
    A scalar specifying the follow-up period of the trial, that is, all patients are followed during this period right after the accrual period.
  }
  \item{method}{
    Character specifying which method would be used to estimate event probabilities. Choices are \code{KM} for the Kaplan-Meier estimator and \code{NNE} for the nearest neighbor estimator described in Heagerty et al (2000). Algorithm of NNE was adapted from the code in the \code{survivalROC} R package. If \code{method = "NNE"} is specified, only fixed-length trials (instead of accrual and follow-up) and \code{cost.keeping} (instead of \code{cost.unit.keeping}) can be specified, and no standard errors will be estimated. Default is \code{KM}.
  }
  \item{lambda}{
    The smoothing parameter used by method \code{NNE} (see documentation of package \code{survivalROC}). Default is 0.05.
  }
  \item{cost.screening}{
    (Optional) the cost of measuring the biomarker for each patient to determine eligibility for the trial.
  }
  \item{cost.keeping}{
    (Optional) the cost of enrolling and retaining a patient in a trial. If specified for a trial with accrual and follow-up periods, it is treated as the "average" cost for one patient under such design, and is used for all patients in the trial.
  }
  \item{cost.unit.keeping}{
    (Optional) the cost of enrolling and retaining a patient per unit time in the trial. Can be specified as an alternative to \code{cost.keeping}, assuming that patients who experienced a clinical event no longer cost in the trial.
  }
  \item{power}{
    The power (probability of rejecting the null hypothesis given that it is false) for sample size calculation.
  }
  \item{alpha}{
    The type I error rate (probability of falsely rejecting the null given that it is true) for sample size calculation.
  }
  \item{one.sided}{
    Logical indicating whether the alternative hypothesis is one-sided (instead of two-sided). Default is \code{FALSE}.
  }
  \item{selected.biomarker.quantiles}{
    Numeric vector specifying the quantiles of the biomarker measured in controls that will be used to screen trial participants. Default is 0, 0.05, ..., 0.95. All entries must be between at least 0 and less than 1.
  }
  \item{do.bootstrap}{
    Logical specifying whether bootstrap standard errors should be calculated for a trial with accrual and follow-up periods. Default is \code{FALSE}. (Standard errors will always be calculated for fixed-length trials.)
  }
  \item{n.bootstrap}{
    Number of bootstrap samples for standard error estimation.
  }
  \item{seed}{
    Specify the random seed used to generate bootstrap samples.
  }
  \item{print.summary.tables}{
    Logical specifying whether a table of clinical trial metrics should be printed.
  }
}
\value{
  Returns a list of the following items:
  \item{summary.table}{
    A matrix containing all clinical trial metrics that were calculated (listed below) at each enrichment level.
  }
  \item{event.prob}{
    Vector(s) of estimated event probabilities at each enrichment level. If more than one vector was presented, each column corresponds to a specified length of trial.
  }
  \item{event.prob.se}{
    Estimated standard errors of event probabilities at each enrichment level for each trial length.
  }
  \item{n.patients}{
    Vector(s) of clinical trial sample size required.
  }
  \item{n.patients.se}{
    Estimated standard errors of the sample sizes.
  }
  \item{n.screened}{
    The number of patients that need to be screened to enroll the trial.
  }
  \item{n.screened.se}{
    Estimated standard errors of \code{n.screened}.
  }
  \item{cost}{
    The estimated total cost of the trial, combining screening costs and the total cost of patients in the trial.
  }
  \item{cost.se}{
    Estimated standard errors of \code{cost}.
  }
  \item{cost.reduction}{
    The reduction in total cost comparing an enriched trial to an unenriched trial. A positive number indicates that an enriched trial would cost less than an unenriched one.
  }
  \item{cost.reduction.se}{
    Estimated standard errors of \code{cost.reduction} from bootstrap.
  }
  \item{response}{
    The response (a Survival object) specified in argument \code{formula}.
  }
  \item{biomarker}{
    The biomarker, or the composite biomarker calculated via Cox regression.
  }
All other quantities are the same as the input by user.
}
\examples{
require(survival)
data(SurvMarkers)

SurvMarkers$surv <- Surv(SurvMarkers$time, SurvMarkers$event)

## Using a composite biomarker (combining x1 and x2) to plan a trial
## with an accrual period of 12 months and a follow-up period of 36 months.
## The cost for keeping a patient in the trial is $300/month.
## We wish to detect a hazard ratio of 0.8 with a two-sided test,
## with type I error 0.05 and power 0.9.
## The Kaplan-Meier method is used to calculate event probabilities.
rslt <- surv_enrichment(formula = surv~x1+x2, data = SurvMarkers, hr = 0.8, a=12, f=36,
                         cost.screening = 300, cost.keeping = NULL, cost.unit.keeping = 300,
                         method = "KM", power = 0.9, alpha = 0.05, one.sided = FALSE,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.9, by = 0.1),
                         do.bootstrap = FALSE, print.summary.tables = FALSE)
}
\references{
  Heagerty, Patrick J., Thomas Lumley, and Margaret S. Pepe (2000) <doi:10.1111/j.0006-341x.2000.00337.x> "Time-dependent ROC curves for censored survival data and a diagnostic marker"

  Heagerty, Patrick J., Paramita Saha-Chaudhuri, and Maintainer Paramita Saha-Chaudhuri (2013) <https://cran.r-project.org/web/packages/survivalROC/survivalROC.pdf> "Package 'survivalROC'."
}
