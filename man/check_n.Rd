% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_n.R
\name{check_n}
\alias{check_n}
\title{Checks if enough data is available}
\usage{
check_n(trend_variable, min_n = 7)
}
\arguments{
\item{trend_variable}{Vector with values}

\item{min_n}{Minimum amount of observations for trenda}
}
\value{
logical
}
\description{
Returns FALSE for less than 7 observations
}
