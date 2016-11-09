#--------------------#
# Plotting functions #
#--------------------#

# Purpose:
# Defines a set of plotting functions to supplement the
# distribution/density/quantile functions.

# Index
# Lookup - 01:  dist_curve

# Lookup - 01
#' Plot Distribution Function Curves
#'
#' Draws a line for a specified distribution function.
#'
#' @param prm a vector of parameters (distribution dependent).
#' @param ch the choice to condition on (i.e. 0 or 1).
#' @param dist a string indicating the distribution function to
#'   plot. Can be
#'   \describe{
#'     \item{\code{wr}}{The Wald race model.}
#'     \item{\code{ig}}{The inverse gaussian distribution.}
#'   }
#' @param b a vector giving the starting and end point of the curve,
#'   and the number of points used to approximate the curve.
#' @param rt an optional vector of response times to plot the function
#'   over.
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'   }
#' @param ... additional plotting parameters.
#' @return Forthcoming
#' @examples
#' # Forthcoming
#' @export
dist_curve = function( prm, ch = 1, dist = 'wr', b = c(0,2,100),
                       rt = NULL, opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = T else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip )

  if ( length(rt) == 0 ) {
    v = seq( b[1], b[2], length = b[3] )
  } else v = sort( rt )

  # Wald race model
  if ( dist == 'wr' ) {

    # Define default parameters
    pv = c( k1 = 1, xi1 = 4, tau1 = 0,
            k0 = 1, xi0 = 1, tau0 = 0,
            s1 = 1, s0 = 1, rl = 0 )
    nms = names( pv )

    # Extract input
    sel = names( prm )
    if ( length( sel ) == 0 ) {
      names( prm ) = nms[1:length(prm)]
    }
    pv[ names(prm) ] = prm;

    # Joint distribution curve
    p = pwaldrace( v, ch, pv[1], pv[2], pv[3],
                   pv[4], pv[5], pv[6], pv[7],
                   pv[8], pv[9] )

    # If conditional instead of joint
    if (!jnt) {
      adj = pwaldrace( Inf, ch, pv[1], pv[2], pv[3],
                       pv[4], pv[5], pv[6], pv[7],
                       pv[8], pv[9] )
      p = p/adj;
    }

  }

  if ( dist == 'ig' ) {

    pv = c( kappa = 1, xi = 4, sigma = 1 )
    nms = names( pv )
    sel = names( prm )
    if ( length( sel ) == 0 ) {
      names( prm ) = nms[1:length(prm)]
    }
    pv[ names(prm) ] = prm;

    p = pinvgauss( v, pv[1], pv[2], pv[3] )

  }

  if (draw) {
    if (flip) lines( v, -p, ... ) else lines( v, p, ... )
  }
}

