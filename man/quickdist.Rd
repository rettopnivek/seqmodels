% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_functions.R
\name{quickdist}
\alias{quickdist}
\title{Create a 'seqmodels' Object for Quick Plotting}
\usage{
quickdist(dist, type = "PDF", prm = NULL, x = NULL, b = NULL)
}
\arguments{
\item{dist}{a character string indicating the type of distribution to
plot.}

\item{type}{the type of function to compute over the set of plotting
values (either the density, distribution, quantile, or hazard
function).}

\item{prm}{a vector of parameters. If null, default values are used.
A named vector can be passed in to control the order in which
parameter inputs are extracted.}

\item{x}{an optional vector of response times or probabilities to
plot over.}

\item{b}{an optional vector giving the lower and upper limits as
well as the number of increments used to generate a sequence of
plotting values.}
}
\value{
An object of class 'seqmodels', a list consisting of...
}
\description{
A function to create a 'seqmodels' object with methods for the
\code{lines} and \code{points} functions, allowing for quick
plotting of density, distribution, quantile, and hazard functions
over a variety of distributions.
}
\examples{
# Forthcoming

}
