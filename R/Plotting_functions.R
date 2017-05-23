#--------------------#
# Plotting functions #
#--------------------#

# Purpose:
# Defines a set of plotting functions to supplement the
# distribution/density/quantile functions.

# Index
# Lookup - 01:  quickdist
# Lookup - 02:  plot.seqmodels
# Lookup - 03:  lines.seqmodels
# Lookup - 04:  points.seqmodels

# Lookup - 01
#' Create a 'seqmodels' Object for Quick Plotting
#'
#' A function to create a 'seqmodels' object with methods for the
#' \code{lines} and \code{points} functions, allowing for quick
#' plotting of density, distribution, quantile, and hazard functions
#' over a variety of distributions.
#'
#' @param dist a character string indicating the type of distribution to
#'   plot.
#' @param type the type of function to compute over the set of plotting
#'   values (either the density, distribution, quantile, or hazard
#'   function).
#' @param prm a vector of parameters. If null, default values are used.
#'   A named vector can be passed in to control the order in which
#'   parameter inputs are extracted.
#' @param x an optional vector of response times or probabilities to
#'   plot over.
#' @param b an optional vector giving the lower and upper limits as
#'   well as the number of increments used to generate a sequence of
#'   plotting values.
#'
#' @return An object of class 'seqmodels', a list consisting of...
#'
#' @examples
#' # Forthcoming
#'
#' @export
quickdist = function( dist, type = 'PDF', prm = NULL,
                      x = NULL, b = NULL ) {

  stop_f = T # Check if the supplied label matches a distribution

  dist = dist_determine( dist )

  # Exponentially modified Gaussian
  if ( dist == 'emg' ) {
    stop_f = F
    out = dist_emg( prm, type, x, b )
  }
  if ( dist == 'sig' ) {
    stop_f = F
    out = dist_sig( prm, type, x, b )
  }
  if ( dist == 'wr' ) {
    stop_f = F
    out = dist_wr( prm, type, x, b )
  }
  if ( dist == 'wp' ) {
    stop_f = F
    out = dist_wp( prm, type, x, b )
  }
  if ( dist == 'n' ) {
    stop_f = F
    out = dist_n( prm, type, x, b )
  }
  if ( dist == 'ga' ) {
    stop_f = F
    out = dist_ga( prm, type, x, b )
  }
  if ( dist == 'we' ) {
    stop_f = F
    out = dist_we( prm, type, x, b )
  }
  if ( dist == 'ln' ) {
    stop_f = F
    out = dist_ln( prm, type, x, b )
  }
  if ( dist == 'b' ) {
    stop_f = F
    out = dist_b( prm, type, x, b )
  }

  if ( stop_f ) {
    stop( 'Please indicate an appropriate RT distribution.',
          call. = FALSE )
  }

  # Create a 'seqmodels' class for plotting
  class( out ) = 'seqmodels'

  return( out )
}

# Lookup - 02
#' Generate a Blank Plot for a 'seqmodels' Object
#'
#' @param object forthcoming.
#' @param weight forthcoming.
#' @param ... forthcoming.
#'
#' @return Forthcoming
#'
#' @export
plot.seqmodels = function( x, weight = NULL, ... ) {

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
#' Draw a Line Based on a 'seqmodels' Object
#'
#' @param object forthcoming.
#' @param ch forthcoming.
#' @param weight forthcoming.
#' @param ... forthcoming.
#'
#' @return Forthcoming
#'
#' @export
lines.seqmodels = function( object, ch = 1, weight = NULL, ... ) {

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
#' Add a Set of Points Based on a 'seqmodels' Object
#'
#' @param object forthcoming.
#' @param ch forthcoming.
#' @param weight forthcoming.
#' @param ... forthcoming.
#'
#' @return Forthcoming
#'
#' @export
points.seqmodels = function( object, ch = 1, weight = NULL, ... ) {

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
