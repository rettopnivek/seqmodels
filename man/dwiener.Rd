% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dwiener}
\alias{dwiener}
\alias{pwiener}
\alias{qwiener}
\alias{rwiener}
\title{Two-boundary Wiener Process for Choice and Response Times}
\usage{
dwiener(
  rt,
  ch,
  alpha,
  theta,
  xi,
  tau,
  sigma = as.numeric(c(1)),
  ln = FALSE,
  joint = TRUE,
  eps = 1e-29,
  parYes = TRUE
)

pwiener(
  rt,
  ch,
  alpha,
  theta,
  xi,
  tau,
  sigma = as.numeric(c(1)),
  ln = FALSE,
  joint = TRUE,
  lower_tail = TRUE,
  eps = 1e-29,
  parYes = TRUE
)

qwiener(
  p,
  ch,
  alpha,
  theta,
  xi,
  tau,
  sigma = as.numeric(c(1)),
  joint = FALSE,
  eps = 1e-29,
  bounds = 3,
  em_stop = 20,
  err = 1e-08,
  parYes = TRUE
)

rwiener(
  n,
  alpha,
  theta,
  xi,
  tau,
  sigma = as.numeric(c(1)),
  eps = 1e-29,
  bounds = 5,
  em_stop = 30,
  err = 1e-16,
  parYes = TRUE
)
}
\arguments{
\item{rt}{a vector of responses times ( \code{rt} > 0 ).}

\item{ch}{a vector of accuracy/choice values ( \code{ch} = {0,1} ).}

\item{alpha}{a vector of upper boundaries at which the evidence
accumulation terminations.}

\item{theta}{a vector of proportions determining the starting
point for the evidence accumulation, where the starting point
\eqn{\zeta} = \code{alpha}*\code{theta} ( 0 \eqn{\ge} \code{theta}
\eqn{\ge} 1 ).}

\item{xi}{a vector of drift rates, the rate of evidence accumulation
( \code{xi} > 0 ).}

\item{tau}{a vector of residual latencies for the non-decision
component ( \code{tau} > 0 ).}

\item{sigma}{a vector giving the coefficients of drift (also known as
within-trial variability; \code{sigma} > 0 ).}

\item{ln}{logical; if \code{TRUE}, probabilities are given as
log(p).}

\item{joint}{logical; if \code{FALSE} the conditional density
(normalized to integrate to one) is returned. Otherwise, the
joint density (integrating to the choice probability) is
returned.}

\item{eps}{the margin of error for the infinite sums being calculated.}

\item{parYes}{logical; if \code{TRUE} the code is run in parallel.}

\item{lower_tail}{logical; if \code{TRUE} (default), probabilities
are \eqn{P(X \le x)} otherwise \eqn{P( X > x)}.}

\item{p}{a vector of probabilities.}

\item{bounds}{upper limit of the quantiles to explore
for the approximation via linear interpolation.}

\item{em_stop}{the maximum number of iterations to attempt to
find the quantile via linear interpolation.}

\item{err}{the number of decimals places to approximate the
cumulative probability during estimation of the quantile function.}

\item{n}{the number of draws for random generation.}
}
\value{
\code{dwiener} gives the density, \code{pwiener} gives the
distribution function, \code{qwiener} approximates the quantile
function, and \code{rwiener} generates random deviates.

The length of the result is determined by \code{n} for \code{rwiener},
and is the maximum of the length of the numerical arguments for
the other functions.

The numerical arguments other than \code{n} are recycled to the
length of the result.
}
\description{
Density, distribution, quantile, and random generation functions
for a two-boundary wiener process that can be applied to
choice and response times (e.g., Luce, 1986; Ratcliff, 1978).
\code{alpha} refers to the boundary separation, \code{theta}
refers to the proportion governing the start point of accumulation,
\code{xi} refers to the rate of evidence accumulation (drift rate),
\code{tau} refers to the residual latency (e.g., motor and
encoding processes), and \code{sigma} refers to the within-trial
variability of evidence accumulation (the coefficient of drift;
typically set to 1 or 0.1).
}
\section{Details}{

The density function is based on the implementation of Navarro
and Fuss (2009). The distribution function is based on the
implementation of Blurton et al. (2012).

A linear interpolation approach is used to approximate the
quantile function and to random deviates by estimating the
inverse of the cumulative distribution function via an
iterative procedure. When the precision of this estimate is
set to 8 decimal places, the approximation will be typically
accurate to about half of a millisecond.
}

\section{References}{


Blurton, S. P., Kesselmeier, M., & Gondan, M. (2012). Fast and
  accurate calculations for cumulative first-passage time distributions
  in Wiener diffusion models. Journal of Mathematical Psychology,
  56, 470-475.

Luce, R. D. (1986). Response times: Their role in inferring
  elementary mental organization. New York, New York: Oxford University
  Press.

Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations
  for first-passage times in Wiener diffusion models. Journal of
  Mathematical Psychology, 53, 222-230.

Ratcliff, R. (1978). A theory of memory retrieval. Psychological
  review, 85, 59 - 108.
}

\examples{
# Density
dwiener( rt = 0.6, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
  xi = 1.0, tau = 0.3 )
# Distribution function
pwiener( rt = 0.6, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
  xi = 1.0, tau = 0.3 )
# Choice probabilities
pwiener( rt = Inf, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
  xi = 1.0, tau = 0.3 )
# Quantile function (Accurate to ~4 decimal places)
round( qwiener( p = .3499, ch = 1, alpha = 1.6, theta = 0.5,
  xi = 1.0, tau = 0.3 ), 4 )
# For quantiles based on joint distribution, re-weight input 'p'
# based on choice probabilities
prob <- pwiener( rt = Inf, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
  xi = 1.0, tau = 0.3 )
round( qwiener( p = .3499 * prob, ch = c( 1, 0), alpha = 1.6,
  theta = 0.5, xi = 1.0, tau = 0.3 ), 4 )
# Simulate values
sim <- rwiener( n = 100, alpha = 0.8, theta = 0.6,
  xi = 0.0, tau = 0.3 )

# Plotting
layout( matrix( 1:4, 2, 2, byrow = T ) )
# Parameters
prm <- c( a = 1.2, z = .4, v = 1.0, t0 = 0.3 )
# Density
obj <- quickdist( 'wp', 'PDF', prm )
plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
# CDF
obj <- quickdist( 'wp', 'CDF', prm )
plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
# Quantiles
obj <- quickdist( 'wp', 'QF', prm, x = seq( .2, .8, .2 ) )
plot( obj ); prb = seq( .2, .8, .2 )
abline( h = prb, lty = 2 )
# Conditional, not joint
lines( obj, type = 'b', pch = 19, weight = 1 )
lines( obj, ch = 0, type = 'b', pch = 17, lty = 2, weight = 1 )
# Hazard function
obj <- quickdist( 'wp', 'HF', prm )
plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )

}
