\name{mphgr_run}
\alias{mphgr_run}
\title{Run microbiome personalized hyper-graph ranking}
\description{
Calls the bundled Python implementation using \pkg{reticulate}.
Input microbiome matrices must already be preprocessed.
}
\usage{
mphgr_run(data_dir, data_type, d, Q, prefix_data, stages = NULL)
}
\arguments{
\item{data_dir}{Folder containing \file{control.txt} and numbered test files.}
\item{data_type}{\code{"1"} for longitudinal samples; \code{"2"} for stage-grouped samples.}
\item{d}{Minimum association count with \code{abs(PCC) >= 0.3}.}
\item{Q}{Differential correlation cutoff for hyper-graph construction.}
\item{prefix_data}{Test file prefix, e.g. \code{"test2_data"} for \file{test2_data1.txt}.}
\item{stages}{Named list of stages and numeric file IDs for mode \code{"2"}. The first stage is the normalization reference.}
}
\details{
Input files must have taxa in rows and samples in columns. The function
writes \file{Score.csv} and \file{d.csv} to \code{data_dir}. Mode \code{"2"}
retains the Python stage-score plot. Python packages \pkg{numpy}, \pkg{pandas},
and \pkg{matplotlib} must be available to \pkg{reticulate}.
}
\value{
Invisibly returns a list: \code{normalized_scores} for mode \code{"1"},
or \code{stage_scores} for mode \code{"2"}.
}
\examples{
\dontrun{
library(mPHGR)
# PTB longitudinal data
mphgr_run("C:/mPHGR/PTB10013", "1", 30, 0.5, "test2_data")
# Gastric cancer stages
mphgr_run("C:/mPHGR/GC", "2", 145, 0.7, "test3_data",
          stages = list(SG = 1:7, AG = 8:17, GIN = 18:22, GC = 23:37))
}
}
