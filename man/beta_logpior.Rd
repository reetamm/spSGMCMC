% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{beta_logpior}
\alias{beta_logpior}
\title{log-prior for regression coefficients \eqn{\beta}
We consider a standard normal prior on all the components}
\usage{
beta_logpior(beta, var = 1)
}
\arguments{
\item{beta, }{the current value of the regression parameters}
}
\value{
a scalar i.e the log-prior at beta
}
\description{
log-prior for regression coefficients \eqn{\beta}
We consider a standard normal prior on all the components
}
