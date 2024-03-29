% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple_rs.R
\name{simple_rs_probabilities}
\alias{simple_rs_probabilities}
\title{Inclusion Probabilities: Simple Random Sampling}
\usage{
simple_rs_probabilities(
  N,
  prob = NULL,
  prob_unit = NULL,
  check_inputs = TRUE,
  simple = TRUE
)
}
\arguments{
\item{N}{The number of units. N must be a positive integer. (required)}

\item{prob}{prob is the probability of being sampled must be a real number between 0 and 1 inclusive, and must be of length 1. (optional)}

\item{prob_unit}{prob is the probability of being sampled must be a real number between 0 and 1 inclusive, and must be of length N. (optional)}

\item{check_inputs}{logical. Defaults to TRUE.}

\item{simple}{logical. internal use only.}
}
\value{
A vector length N indicating the probability of being sampled.
}
\description{
Inclusion Probabilities: Simple Random Sampling
}
\examples{
probs <- simple_ra_probabilities(N = 100)
table(probs)

probs <- simple_ra_probabilities(N = 100, prob = 0.3)
table(probs)

}
