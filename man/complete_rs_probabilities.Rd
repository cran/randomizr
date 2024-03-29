% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complete_rs.R
\name{complete_rs_probabilities}
\alias{complete_rs_probabilities}
\title{Inclusion Probabilities: Complete Random Sampling}
\usage{
complete_rs_probabilities(
  N,
  n = NULL,
  n_unit = NULL,
  prob = NULL,
  prob_unit = NULL,
  check_inputs = TRUE
)
}
\arguments{
\item{N}{The number of units. N must be a positive integer. (required)}

\item{n}{Use for a design in which exactly n units are sampled. (optional)}

\item{n_unit}{unique(n_unit) will be passed to \code{n}. Must be the same for all units (optional)}

\item{prob}{Use for a design in which either floor(N*prob) or ceiling(N*prob) units are sampled. The probability of being sampled is exactly prob because with probability 1-prob, floor(N*prob) units will be sampled and with probability prob, ceiling(N*prob) units will be sampled. prob must be a real number between 0 and 1 inclusive. (optional)}

\item{prob_unit}{unique(prob_unit) will be passed to the prob argument and must be the same for all units.}

\item{check_inputs}{logical. Defaults to TRUE.}
}
\value{
A vector length N indicating the probability of being sampled.
}
\description{
Inclusion Probabilities: Complete Random Sampling
}
\examples{
probs <- complete_rs_probabilities(N = 100)
table(probs)

probs <- complete_rs_probabilities(N = 100, n = 50)
table(probs)

probs <- complete_rs_probabilities(N=100, prob = .3)
table(probs)

}
