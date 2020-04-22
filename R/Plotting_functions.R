#--------------------#
# Plotting functions #
#--------------------#

# Purpose:
# Defines a set of plotting functions to supplement the
# distribution/density/quantile functions.

# Index
# Lookup - 01:  quickdist
# Lookup - 02:  plot.quickdist
# Lookup - 03:  lines.quickdist
# Lookup - 04:  points.quickdist
# Lookup - 05:  is.quickdist

# Lookup - 01
#' Create a 'quickdist' Object for Quick Plotting
#'
#' A function to create a \code{quickdist} object with methods for the
#' \code{plot}, \code{lines}, and \code{points} functions,
#' allowing for quick plotting of density, distribution, quantile,
#' and hazard functions over a variety of distributions.
#'
#' @param dist a character string indicating the type of distribution to
#'   plot.
#' @param type the type of function to compute over the set of plotting
#'   values, either the density ('PDF', 'pdf', 'density', 'Density',
#'   'd'), distribution ('CDF', 'cdf', 'Distribution', 'distribution',
#'   'df', 'DF'), quantile ('QF', 'qf', 'QPE', 'qpe',
#'   'quantile', 'quantiles', 'Quantile', 'Quantiles', 'q'), or hazard
#'   function ('HF', 'hf', 'h', 'Hazard', 'hazard').
#' @param prm a vector of parameters. If null, default values are used.
#'   A named vector can be passed in to control the order in which
#'   parameter inputs are extracted.
#' @param x an optional vector of response times or probabilities to
#'   plot over.
#' @param b an optional vector giving the lower and upper limits as
#'   well as the number of increments used to generate a sequence of
#'   plotting values.
#'
#' @details Currently supported distributions include...
#' \itemize{
#'   \item Exponentially modified Gaussian ('emg', 'EMG');
#'   \item Shifted inverse Gaussian, or Wald
#'     ('sig', 'SIG', 'wald', 'Wald', 'ig', 'IG' );
#'   \item Two-boundary Wiener process, or Wiener
#'     diffusion ('wp', 'WP', 'wd', 'WD', 'wiener', 'Wiener');
#'   \item Gaussian, or Normal ('n', 'N', 'gauss', 'Gauss',
#'     'normal', 'Normal', 'gaussian', 'Gaussian' );
#'   \item Gamma ('ga', 'gamma', 'Gamma');
#'   \item Weibull ('wei', 'weibull', 'Weibull');
#'   \item Log-normal ('ln', 'LN', 'log-normal', 'Log-normal');
#'   \item Beta ('b', 'B', 'beta', 'Beta');
#'   \item Exponential ('exp', 'Exp', 'exponential',
#'     'Exponential');
#' }
#'
#' Once a \code{quickdist} object has been created, a blank plot
#' can be generated via the \code{plot} command. The figure can
#' then be updated via the \code{lines} or \code{points} functions,
#' making it easy to add multiple functions to a single figure.
#'
#' @return An object of class 'quickdist', a list consisting of...
#' \itemize{
#'   \item \code{d_pv} - x and y-axis values for the density
#'     function;
#'   \item \code{p_pv} - x and y-axis values for the distribution
#'     function;
#'   \item \code{q_pv} - x and y-axis values for the quantile
#'     function;
#'   \item \code{h_pv} - x and y-axis values for the hazard
#'     function;
#'   \item \code{type} - a numeric code indicating the
#'     type of function to plot.
#'   \item \code{prm} - a named vector with the parameter
#'     values for the specified distribution.
#'   \item \code{dist} - A character string with the
#'     abbreviated label for the specified distribution.
#'   \item \code{ch} - A numeric code indicating the
#'     choice for the conditional distribution (for
#'     two-choice models).
#' }
#'
#' @examples
#' # Plot density function for shifted inverse Gaussian
#' obj <- quickdist( 'SIG' )
#' plot( obj ); lines( obj )
#' # Add line for density function with new values
#' # for kappa, xi, and tau (threshold, drift, and
#' # shift)
#' obj_2 <- quickdist( 'SIG', prm = c( .5, 2, .3 ) )
#' lines( obj_2, col = 'red', lty = 2 )
#'
#' # Plot distribution function for normal distribution
#' # over custom range and number of points used to
#' # approximate line
#' obj <- quickdist( 'Normal', type = 'CDF', b = c(-4,4,10) )
#' plot( obj ); lines( obj ); points( obj, pch = 19 )
#'
#' # Plot predicted versus observed quantile functions
#'
#' # Simulate data for create observed function
#' prb <- seq( .1, .9, .1 ) # 10% to 90% in 10% increments
#' # Simulate data from Beta( 5, 2 )
#' x <- rbeta( 1000, 5, 2 )
#' q_obs <- list( x = quantile( x, prob = prb ) )
#' q_obs$y <- prb
#'
#' # Examine quantile functions for two different predictions
#' # Predicted under Beta(1,1)
#' prd_1 <- quickdist( 'b', type = 'QF', x = prb, prm = c( 1, 1 ) )
#' # Predicted under Beta(5,2)
#' prd_2 <- quickdist( 'Beta', type = 'QF', x = prb, prm = c( 5, 2 ) )
#'
#' plot( prd_1 ); points( prd_1, col = 'blue', pch = 19 )
#' points( prd_2, col = 'red', pch = 15 )
#' # Observed
#' points( q_obs$x, q_obs$y, cex = 1.5 )
#'
#' # Plot hazard function for exponential distribution
#' # (demonstrates memoryless property)
#' obj <- quickdist( 'Exponential', type = 'HF' )
#' plot( obj ); lines( obj )
#' # Contrast with gamma distribution
#' obj <- quickdist( 'Gamma', type = 'HF', prm = c( 1.1, 1) )
#' lines( obj, lty = 2 )
#'
#' @export

quickdist = function( dist, type = 'PDF', prm = NULL,
                      x = NULL, b = NULL ) {

  stop_f = T # Check if the supplied label matches a distribution

  dist = dist_determine( dist )

  # Exponentially modified Gaussian
  if ( dist %in% c( 'emg', 'EMG' ) ) {
    stop_f = F
    out = dist_emg( prm, type, x, b )
  }
  # Shifted inverse Gaussian
  if ( dist %in% c( 'sig', 'SIG', 'wald', 'Wald', 'ig', 'IG' ) ) {
    stop_f = F
    out = dist_sig( prm, type, x, b )
  }
  #if ( dist == 'wr' ) {
  #  stop_f = F
  #  out = dist_wr( prm, type, x, b )
  #}
  # Two-boundary Wiener process
  if ( dist %in% c( 'wp', 'WP', 'wd', 'WD', 'wiener', 'Wiener' ) ) {
    stop_f = F
    out = dist_wp( prm, type, x, b )
  }
  # Gaussian
  if ( dist %in% c( 'n', 'N', 'gauss', 'Gauss', 'gaussian', 'Gaussian',
                    'normal', 'Normal', 'norm', 'Norm' ) ) {
    stop_f = F
    out = dist_n( prm, type, x, b )
  }
  # Gamma
  if ( dist %in% c( 'ga', 'gamma', 'Gamma' ) ) {
    stop_f = F
    out = dist_ga( prm, type, x, b )
  }
  # Weibull
  if ( dist %in% c( 'wei', 'weibull', 'Weibull' ) ) {
    stop_f = F
    out = dist_we( prm, type, x, b )
  }
  # Log-normal
  if ( dist %in% c( 'ln', 'LN', 'log-normal', 'Log-normal' ) ) {
    stop_f = F
    out = dist_ln( prm, type, x, b )
  }
  # Beta
  if ( dist %in% c( 'b', 'B', 'beta', 'Beta' ) ) {
    stop_f = F
    out = dist_b( prm, type, x, b )
  }
  # Exponential
  if ( dist %in% c( 'exp', 'exp', 'exponential', 'Exponential' ) ) {
    stop_f = F
    out = dist_exp( prm, type, x, b )
  }

  if ( stop_f ) {
    stop( paste0(
      'Please indicate an appropriate distribution. Try ?quickdist ',
      'to see supported distributions.' ),
          call. = FALSE )
  }

  # Create a 'quickdist' class for plotting
  class( out ) = 'quickdist'

  return( out )
}

# Lookup - 02
#' @rdname quickdist
#' @export

plot.quickdist = function( x, weight = NULL, ... ) {

  # Extract 'seqmodels' object
  object = x

  # Default labels
  labels = c( 'Time', 'Density function' )
  if ( object$dist == 'b' ) labels[1] = 'Probability'
  if ( object$dist == 'n' ) labels[1] = 'Z-scores'

  # Extract x and y values
  if ( object$ch == 1 ) {

    if ( is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      xv = object$d_pv$x
      yv = object$d_pv$y
    }
    if ( object$type == 2 ) {
      xv = object$p_pv$x
      yv = object$p_pv$y
      labels[2] = 'Distribution function'
    }
    if ( object$type == 3 ) {
      xv = object$q_pv$x
      yv = object$q_pv$y
      labels[2] = 'Quantile function'
    }
    if ( object$type == 4 ) {
      xv = object$h_pv$x
      yv = object$h_pv$y
      labels[2] = 'Hazard function'
    }
  }

  if ( object$ch == 2 ) {

    if ( object$type != 3 && is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      xv = object$d_pv$x
      cmp = sapply( object$d_pv[-1], function(x) max( unlist( x ) ) )
      if ( max( cmp ) == cmp[1] ) yv = object$d_pv$y1 else
        yv = object$d_pv$y0
    }
    if ( object$type == 2 ) {
      xv = object$p_pv$x
      cmp = sapply( object$p_pv[-1], function(x) max( unlist( x ) ) )
      if ( max( cmp ) == cmp[1] ) yv = object$p_pv$y1 else
        yv = object$p_pv$y0
      labels[2] = 'Distribution function'
    }
    if ( object$type == 3 ) {
      if ( is.null( weight ) ) {
        weight = max( object$p1, 1 - object$p1 )
      }
      yv = object$q_pv$y
      cmp = sapply( object$q_pv[-1], function(x) max( unlist( x ) ) )
      if ( max( cmp ) == cmp[1] ) xv = object$q_pv$x1 else
        xv = object$q_pv$x0
      labels[2] = 'Quantile function'
    }
    if ( object$type == 4 ) {
      xv = object$h_pv$x
      cmp = sapply( object$h_pv[-1], function(x) max( unlist( x ) ) )
      if ( max( cmp ) == cmp[1] ) yv = object$h_pv$y1 else
        yv = object$h_pv$y0
      labels[2] = 'Hazard function'
    }
  }

  # Convert variables to list
  arguments = list( ... )
  int = arguments[['int']]
  if ( is.null( int ) ) int = .5

  # Extract names of arguments
  arg_names = names( arguments )
  rm( arguments )

  if ( any( 'type' %in% arg_names ) ) {
    stop( "Please exclude variable 'type'",
          call. = FALSE )
  }

  # Lower and upper limits based on
  # 'seqmodels' object
  xl = lower_upper( int, xv )
  yl = lower_upper( int, yv )
  # x and y-axis labels
  xlb = labels[1];
  ylb= labels[2];

  # Check for the variables
  # 'xlab', 'ylab'

  # No xlab, ylab
  if ( !any( 'xlab' %in% arg_names ) &
       !any( 'ylab' %in% arg_names ) ) {
    plot( x = xl, y = yl, type = 'n',
          xlab = xlb, ylab = ylb, ... )
  }

  # No xlab
  if ( !any( 'xlab' %in% arg_names ) &
       any( 'ylab' %in% arg_names ) ) {
    plot( x = xl, y = yl, type = 'n',
          ylab = ylb, ... )
  }

  # No ylab
  if ( any( 'xlab' %in% arg_names ) &
       !any( 'ylab' %in% arg_names ) ) {
    plot( x = xl, y = yl, type = 'n',
          ylab = ylb, ... )
  }

  # All present
  if ( any( 'xlab' %in% arg_names ) &
       any( 'ylab' %in% arg_names ) ) {
    plot( x = xl, y = yl,
          ylab = ylb, type = 'n', ... )
  }

}

# Lookup - 03
#' @rdname quickdist
#' @export

lines.quickdist = function( object, ch = 1, weight = NULL, ... ) {

  # Extract x and y values
  if ( object$ch == 1 ) {

    if ( is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      x = object$d_pv$x
      y = object$d_pv$y
    }
    if ( object$type == 2 ) {
      x = object$p_pv$x
      y = object$p_pv$y
    }
    if ( object$type == 3 ) {
      x = object$q_pv$x
      y = object$q_pv$y
    }
    if ( object$type == 4 ) {
      x = object$h_pv$x
      y = object$h_pv$y
    }
  }

  if ( object$ch == 2 ) {

    if ( object$type != 3 && is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      x = object$d_pv$x
      if ( ch == 1 ) y = object$d_pv$y1 else y = object$d_pv$y0
    }
    if ( object$type == 2 ) {
      x = object$p_pv$x
      if ( ch == 1 ) y = object$p_pv$y1 else y = object$p_pv$y0
    }
    if ( object$type == 3 ) {
      if ( is.null( weight ) ) {
        if ( ch == 1 ) weight = object$p1
        if ( ch == 0 ) weight = 1 - object$p1
      }
      y = object$q_pv$y
      if ( ch == 1 ) x = object$q_pv$x1 else x = object$q_pv$x0
    }
    if ( object$type == 4 ) {
      x = object$h_pv$x
      if ( ch == 1 ) y = object$h_pv$y1 else y = object$h_pv$y0
    }
  }

  # Draw line
  lines( x, weight * y, ... )
}

# Lookup - 04
#' @rdname quickdist
#' @export

points.quickdist = function( object, ch = 1, weight = NULL, ... ) {

  # Extract x and y values
  if ( object$ch == 1 ) {

    if ( is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      x = object$d_pv$x
      y = object$d_pv$y
    }
    if ( object$type == 2 ) {
      x = object$p_pv$x
      y = object$p_pv$y
    }
    if ( object$type == 3 ) {
      x = object$q_pv$x
      y = object$q_pv$y
    }
    if ( object$type == 4 ) {
      x = object$h_pv$x
      y = object$h_pv$y
    }
  }

  if ( object$ch == 2 ) {

    if ( object$type != 3 && is.null( weight ) ) weight = 1

    if ( object$type == 1 ) {
      x = object$d_pv$x
      if ( ch == 1 ) y = object$d_pv$y1 else y = object$d_pv$y0
    }
    if ( object$type == 2 ) {
      x = object$p_pv$x
      if ( ch == 1 ) y = object$p_pv$y1 else y = object$p_pv$y0
    }
    if ( object$type == 3 ) {
      if ( is.null( weight ) ) {
        if ( ch == 1 ) weight = object$p1
        if ( ch == 0 ) weight = 1 - object$p1
      }
      y = object$q_pv$y
      if ( ch == 1 ) x = object$q_pv$x1 else x = object$q_pv$x0
    }
    if ( object$type == 4 ) {
      x = object$h_pv$x
      if ( ch == 1 ) y = object$h_pv$y1 else y = object$h_pv$y0
    }
  }

  # Draw line
  points( x, weight * y, ... )
}

# Lookup - 05
#' @rdname quickdist
#' @export

is.quickdist <- function(x) inherits(x, "quickdist")
