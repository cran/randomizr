% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/declare_rs.R
\name{declare_rs}
\alias{declare_rs}
\title{Declare a random sampling procedure.}
\usage{
declare_rs(
  N = NULL,
  strata = NULL,
  clusters = NULL,
  n = NULL,
  n_unit = NULL,
  prob = NULL,
  prob_unit = NULL,
  strata_n = NULL,
  strata_prob = NULL,
  simple = FALSE,
  check_inputs = TRUE
)
}
\arguments{
\item{N}{The number of units. N must be a positive integer. (required)}

\item{strata}{A vector of length N that indicates which stratum each unit belongs to.}

\item{clusters}{A vector of length N that indicates which cluster each unit belongs to.}

\item{n}{Use for a design in which n units (or clusters) are sampled. In a stratified design, exactly n units in each stratum will be sampled. (optional)}

\item{n_unit}{Under complete random sampling, must be constant across units. Under stratified random sampling, must be constant within strata.}

\item{prob}{Use for a design in which either floor(N*prob) or ceiling(N*prob) units (or clusters) are sampled. The probability of being sampled is exactly prob because with probability 1-prob, floor(N*prob) units (or clusters) will be sampled and with probability prob, ceiling(N*prob) units (or clusters) will be sampled. prob must be a real number between 0 and 1 inclusive. (optional)}

\item{prob_unit}{Must of be of length N. Under simple random sampling, can be different for each unit or cluster.  Under complete random sampling, must be constant across units. Under stratified random sampling, must be constant within strata.}

\item{strata_n}{Use for a design in which strata_n describes the number of units to sample within each stratum.}

\item{strata_prob}{Use for a design in which strata_prob describes the probability of being sampled within each stratum. Differs from prob in that the probability of being sampled can vary across strata.}

\item{simple}{logical, defaults to FALSE. If TRUE, simple random sampling is used. When \code{simple = TRUE}, please do not specify n or strata_n. When \code{simple = TRUE}, \code{prob} may vary by unit.}

\item{check_inputs}{logical. Defaults to TRUE.}
}
\value{
A list of class "declaration".  The list has five entries:
  $rs_function, a function that generates random samplings according to the declaration.
  $rs_type, a string indicating the type of random sampling used
  $probabilities_vector, A vector length N indicating the probability of being sampled.
  $strata, the stratification variable.
  $clusters, the clustering variable.
}
\description{
Declare a random sampling procedure.
}
\examples{
# The declare_rs function is used in three ways:

# 1. To obtain some basic facts about a sampling procedure:
declaration <- declare_rs(N = 100, n = 30)
declaration

# 2. To draw a random sample:

S <- draw_rs(declaration)
table(S)

# 3. To obtain inclusion probabilities

probs <- obtain_inclusion_probabilities(declaration)
table(probs, S)

# Simple Random Sampling Declarations

declare_rs(N = 100, simple = TRUE)
declare_rs(N = 100, prob = .4, simple = TRUE)

# Complete Random Sampling Declarations

declare_rs(N = 100)
declare_rs(N = 100, n = 30)

# Stratified Random Sampling Declarations

strata <- rep(c("A", "B","C"), times=c(50, 100, 200))
declare_rs(strata = strata)
declare_rs(strata = strata, prob = .5)


# Cluster Random Sampling Declarations

clusters <- rep(letters, times = 1:26)
declare_rs(clusters = clusters)
declare_rs(clusters = clusters, n = 10)

# Stratified and Clustered Random Sampling Declarations

clusters <- rep(letters, times = 1:26)
strata <- rep(NA, length(clusters))
strata[clusters \%in\% letters[1:5]] <- "stratum_1"
strata[clusters \%in\% letters[6:10]] <- "stratum_2"
strata[clusters \%in\% letters[11:15]] <- "stratum_3"
strata[clusters \%in\% letters[16:20]] <- "stratum_4"
strata[clusters \%in\% letters[21:26]] <- "stratum_5"

table(strata, clusters)

declare_rs(clusters = clusters, strata = strata)
declare_rs(clusters = clusters, strata = strata, prob = .3)

}
