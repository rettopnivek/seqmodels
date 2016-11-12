#--------------------#
# Plotting functions #
#--------------------#

# Purpose:
# Defines a set of plotting functions to supplement the
# distribution/density/quantile functions.

# Index
# Lookup - 01:  dist_curve

# Lookup - 01
#' Plot Functions of Distributions
#'
#' Adds a lines/set of points to an existing plot for a
#' desired function (e.g. cumulative density, probability density,
#' hazard) of a response time distribution.
#'
#' @param prm a vector of distribution-dependent parameters.
#' @param rt an optional vector of response times to plot the function
#'   over.
#' @param ch the choice to condition on (i.e. 0 or 1).
#' @param dist a string indicating the distribution function to
#'   plot. Can be
#'   \describe{
#'     \item{\code{wr}}{The Wald race model.}
#'     \item{\code{ig}}{The inverse gaussian distribution.}
#'     \item{\code{emg}}{The exponentially modified gaussian
#'       distribution.}
#'   }
#' @param ver the type of function for the distribution. Can be
#'   \describe{
#'     \item{\code{CDF}}{The cumulative density function.}
#'     \item{\code{PDF}}{The probability density function.}
#'     \item{\code{QPE}}{The quantile probability points.}
#'     \item{\code{HF}}{The hazard function.}
#'   }
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{b}}{The lower and upper boundaries for
#'       the sequence of times and the number of points to
#'       use to approximate the curve.}
#'     \item{\code{prb}}{The quantile probabilities to
#'       calculate.}
#'     \item{\code{pts}}{If false, a line is drawn instead of
#'       a set of points.}
#'     \item{\code{draw}}{If true, the line or set of points is
#'       drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve or set of points
#'       is flipped about the x-axis.}
#'   }
#' @param ... additional plotting parameters.
#' @return A list with the x and y-axis plotting values.
#' @examples
#' # Forthcoming
#' @export

dist_plot = function( prm = NULL, rt = NULL, ch = 1,
                      dist = 'wr', ver = 'CDF',
                      opt = list( ),
                      ... ) {

  # Set options
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$b ) == 0 ) b = c(0,2,100) else b = opt$b
  if ( length( opt$prb ) == 0 ) prb = seq(.1,.9,.2) else
    prb = opt$prb
  if ( length( opt$pts ) == 0 ) pts = F else
    pts = opt$pts

  # If no specified set of response times,
  # create sequence to plot over
  if ( length( rt ) == 0 ) {
    t = seq( b[1], b[2], length = b[3] )
  } else {
    t = rt
  }

  output = dist_calc(t,ch,prm,dist,prb,ver)
  t = output$t; y = output$y;

  if (draw & length(t) != 0 ) {
    if (pts) {
      if (!flip) points( t, y, ... ) else points( t, -y, ... )
    } else {
      if (!flip) lines( t, y, ... ) else lines( t, -y, ... )
    }
  }

  names( output ) = c('x','y')
  if (out) return( output )
}
