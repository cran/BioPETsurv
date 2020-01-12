\name{surv_plot_enrichment}
\alias{surv_plot_enrichment}
\title{Plotting Clinical Trial Metrics for Prognostic Enrichment}
\description{
This function plots summaries of prognostic enrichment of clinical trials with survival outcomes, based on clinical trial metrics estimated by \code{surv_enrichment}.
}
\usage{
surv_plot_enrichment(x, km.quantiles = c(0,0.25,0.5,0.75),
                     km.range = NULL, alt.color = NULL)
}
\arguments{
  \item{x}{
    Object returned by \code{surv_enrichment}.
  }
  \item{km.quantiles}{
    Enrichment levels on which Kaplan-Meier survival estimates (Plot 1) are plotted. Defaults to four quartiles.
  }
  \item{km.range}{
    (Optional) a scalar specifying the range of time for which Kaplan-Meier survival estimates (Plot 1) are plotted. Defaults to the last time point of observation.
  }
  \item{alt.color}{
    (Optional) allows the user to specify the color of curves for clinical trial metrics (Plots 2-6). The length should match the number of trial lengths considered. Defaults to \code{ggplot2} color palette.
  }
}
\value{
A grid containing either the first 4 or 6 plots described below.
  \item{km.plot}{
    The Kaplan-Meier survival curves for specified enrichment levels. The vertical reference line(s) correspond to \code{end.of.trial} or \code{a,f}. This will be presented even if \code{method = "NNE"} was specified.
  }
  \item{prob.plot}{
    The estimated event probability (and 95\% confidence intervals) at each enrichment level.
  }
  \item{ss.plot}{
    The estimated sample size (and confidence intervals) at each enrichment level.
  }
  \item{screen.plot}{
    The estimated number of patients that need to be screened (and confidence intervals) to enroll the trial.
  }
  \item{cost.plot}{
    The estimated total cost of the trial (and confidence intervals).
  }
  \item{reduction.cost.plot}{
    The percentage of reduction in total cost comparing an enriched versus unenriched trial.
  }
  \item{summary}{
    A grid of the first 4 or all 6 plots combined together.
  }
}
\examples{
## Following the example of 'surv_enrichment':
data(SurvMarkers)
library(survival)

SurvMarkers$surv <- Surv(SurvMarkers$time, SurvMarkers$event)
rslt <- surv_enrichment(formula = surv~x1+x2, data = SurvMarkers, hr = 0.8, a=12, f=36,
                         cost.screening = 300, cost.keeping = NULL, cost.unit.keeping = 300,
                         method = "KM", power = 0.9, alpha = 0.05, one.sided = FALSE,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.9, by = 0.1),
                         do.bootstrap = FALSE, print.summary.tables = FALSE)

plots <- surv_plot_enrichment(rslt, km.quantiles = c(0,0.25,0.5,0.75))
}
