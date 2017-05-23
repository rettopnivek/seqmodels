#-------------------#
# Utility functions #
#-------------------#

# Index
# Lookup - 01:  prm_create
# Lookup - 02:  check_version
# Lookup - 03:  create_bounds
# Lookup - 04:  create_x
# Lookup - 05:  dist_determine
# Lookup - 06:  lower_upper

# Lookup - 01
prm_create = function( dv, prm ) {
  # Purpose:
  # A function that takes a set of default values and replaces
  # them with user-supplied values.
  # Arguments:
  # dv  - A named vector of default parameter inputs
  # prm - A user-supplied vector (optionally named) of parameter
  #       inputs
  # Returns:
  # A new vector with select default values swapped out for the
  # user-supplied values.

  # Initialize output to default values
  out = dv

  # Determine number of default parameters
  L = length(dv)

  # If no parameters were provided
  if ( !is.null( prm ) ) {

    # Determine number of inputs
    Ni = length( prm )

    # Check if length of input is too long
    if ( Ni > L ) {
      warn = paste( 'Too many inputs.', '\n',
                    'Only the first', L,
                    'elements will be used.', '\n' )
      warning( warn, call. = FALSE )
      prm = prm[ 1:L ]
      Ni = L
    }

    # Check if names were provided for parameters
    nms = names( prm )
    if ( !is.null( nms ) ) {

      # Check if given names match required names
      chk = nms %in% names( dv )

      if ( !any( chk ) ) {
        # If no names match, return an error
        err = paste(
          'Names of inputs do not match  parameter names',
        '\n', 'Parameter names are:', '\n',
        paste( names( dv ), collapse = ', ' ) )
        stop( err, call. = FALSE )

      } else {
        # If some/all names match

        # If only some names match, return a warning
        if ( !all( chk ) ) {
          warn = paste(
            'Some names of inputs do not match  parameter names',
            '\n', 'Parameter names are:', '\n',
            paste( names( dv ), collapse = ', ' ) )
          warning( warn, call. = FALSE )
        }

        # Replace default parameters with matching named inputs
        out[ nms[chk] ] = prm[ nms[chk] ]
      }

    } else {
      # If no names are provided, replace default values in sequential
      # order

      tmp = numeric( L )
      tmp[ 1:Ni ] = prm
      if ( L > Ni ) tmp[ (Ni+1):L ] = dv[ (Ni+1):L ]
      names( tmp ) = names( out )
      out = tmp
    }

  }

  return( out )
}

# Lookup - 02
check_version = function( type ) {
  # Purpose:
  # A function that compares the input for the type of
  # function (density, distribution, quantile, or
  # hazard ) against a set of common labels to identify
  # the preferred function to evaluate.
  # Arguments:
  # type - A character string
  # Returns:
  # A numerical index indicating the preferred function to
  # evaluate, where...
  # 1 - Density;
  # 2 - Distribution function;
  # 3 - Quantile function;
  # 4 - Hazard function.

  # Initialize output
  out = NULL

  # Compare input against common labels
  if ( is.character( type ) ) {

    for ( i in 1:4 ) {

      if ( i == 1 ) {
        # Check for PDF
        common_names = c( 'PDF', 'pdf', 'density', 'Density',
                          'd' )
      }
      if ( i == 2 ) {
        # Check for CDF
        common_names = c( 'CDF', 'cdf', 'Distribution', 'distribution',
                          'df', 'DF' )
      }
      if ( i == 3 ) {
        # Check for quantile function
        common_names = c( 'QF', 'qf', 'QPE', 'qpe',
                          'quantile', 'quantiles',
                          'Quantile', 'Quantiles',
                          'q' )
      }
      if ( i == 4 ) {
        # Check for hazard function
        common_names = c( 'HF', 'hf', 'h',
                          'Hazard', 'hazard' )
      }

      # Create numerical index
      chk = type %in% common_names
      if ( all( chk ) ) {
        out = i
        break;
      }

    }

    # If input does not match any labels, return an error
    if ( is.null( out ) ) {
      err = paste(
        'Please specify type of function to return.', '\n',
        '(e.g., PDF, CDF, QF, or HF)', '\n' )
      stop( err, call. = FALSE )
    }
  }

  return( out )
}

# Lookup - 03
create_bounds = function( b, ver ) {
  # Purpose:
  # A function that generates the lower and upper limits
  # to plot over and the number of observations/increments
  # to use.
  # Arguments:
  # b   - A vector of up to 3 values, where...
  #       [1] = lower limit
  #       [2] = upper limit
  #       [3] = number of increments
  # ver - The type of function being evaluated
  # Returns:
  # A vector with 3 values for the lower and upper limits and
  # the number of increments, to be passed into the 'seq'
  # function.

  # If no boundary values were passed in
  if ( is.null(b) ) {

    # Density, distribution, and hazard functions
    if ( ver == 1 | ver == 2 | ver == 4 ) {
      b = c( 0, 2, 1000 )
    }

    # Quantile function
    if ( ver == 3 ) {
      b = c( .1, .9, 5 )
    }

    return( b )
  }

  # If a single value is provided, assume it denotes the
  # upper limit (or a single cumulative probability to evaluate)
  if ( length(b) == 1 ) {

    # Density, distribution, and hazard functions
    if ( ver == 1 | ver == 2 | ver == 4 ) {
      if ( b[1] > 0.0 ) b = c( 0, b[1], 1000 ) else {
        stop( 'Upper limit must be greater than 0',
              call. = FALSE )
      }
    }

    # Quantile function
    if ( ver == 3 ) {

      if ( b[1] >= 0.0 & b[1] <= 1.0 )
        b = c( b[1], b[1], 1 ) else {
          stop( 'Please provide a probability', call. = FALSE )
        }
    }

    return( b )
  }

  # If two values are provided, assume it denotes the lower
  # and upper limit
  if ( length(b) == 2 ) {

    # Density, distribution, and hazard functions
    if ( ver == 1 | ver == 2 | ver == 4 ) {
      if ( b[2] > b[1] ) b = c( b[1], b[2], 1000 ) else {
        stop( 'Upper limit must be higher than lower limit',
              call. = FALSE )
      }
    }

    # Quantile function
    if ( ver == 3 ) {
      if ( b[1] >= 0.0 & b[1] <= 1.0 &
           b[2] >= 0.0 & b[2] <= 1.0 &
           b[2] > b[1] ) {
        inc = round( ( b[2] - b[1] )/.1 )
        b = c( b[1], b[2], inc )
      } else {
        stop( 'Please provide two sequential probabilities',
              call. = FALSE )
      }
    }

    return( b )
  }

  # If three values are provided, assume it denotes the lower
  # and upper limit and the number of increments
  if ( length(b) == 3 ) {

    # Density, distribution, and hazard functions
    if ( ver == 1 | ver == 2 | ver == 4 ) {
      if ( b[2] > b[1] ) b = b else {
        stop( 'Upper limit must be higher than lower limit',
              call. = FALSE )
      }
    }

    # Quantile function
    if ( ver == 3 ) {
      if ( b[1] >= 0.0 & b[1] <= 1.0 &
           b[2] >= 0.0 & b[2] <= 1.0 &
           b[2] > b[1] ) {
        inc = round( ( b[2] - b[1] )/.1 )
        b = c( b[1], b[2], inc )
      } else {
        stop( 'Please provide two sequential probabilities',
              call. = FALSE )
      }
    }

    return( b )
  }

  # If more than three values are provided, truncate and
  # warn the user
  if ( length(b) > 3 ) {

    # Truncate
    b = b[1:3]

    warn = paste(
      'Too many inputs governing values for plotting function.',
      '\n', 'Provide up to three values for lower and upper',
      '\n', 'limits and number of observations to plot.', '\n' )
    warning( warn, call. = FALSE )

    # Density, distribution, and hazard functions
    if ( ver == 1 | ver == 2 | ver == 4 ) {
      if ( b[2] > b[1] ) b = b else {
        stop( 'Upper limit must be higher than lower limit',
              call. = FALSE )
      }
    }

    # Quantile function
    if ( ver == 3 ) {
      if ( b[1] >= 0.0 & b[1] <= 1.0 &
           b[2] >= 0.0 & b[2] <= 1.0 &
           b[2] > b[1] ) {
        inc = round( ( b[2] - b[1] )/.1 )
        b = c( b[1], b[2], inc )
      } else {
        stop( 'Please provide two sequential probabilities',
              call. = FALSE )
      }
    }

    return( b )
  }

}

# Lookup - 04
create_x = function( x, b, ver ) {
  # Purpose:
  # A function that generates a set of values to plot
  # over.
  # Arguments:
  # x   - An optional set of user-supplied values
  # b   - A vector giving the lower and upper limits and the
  #       the number of increments for the plotting values
  # ver - The type of function (density, distribution, quantile,
  #       or hazard) to plot
  # Returns:
  # A sequence of plotting values.

  # If set of values is provided
  if ( !is.null( x ) ) {

    # Check that values are appropriate for the given function type
    if ( ver == 3 ) {
      if ( !all( x >= 0.0 ) & !all( x <= 1.0 ) )
        stop( 'Please provide probabilities', call. = FALSE )
    }

    return( x )
  } else {
    b = create_bounds( b, ver )
    x = seq( b[1], b[2], length = b[3] )
    return( x )
  }

}

# Lookup - 05
dist_determine = function( dist ) {
  # Purpose:
  # A function that attempts to determine the distribution
  # that matches the inputted character string.
  # Arguments:
  # dist - An inputted character string.
  # Returns:
  # The internal label the plotting function uses to match
  # distributions.

  out = NULL

  available_dist = matrix( ' ', 9, 2 );
  available_dist[,1] = c(
    'Ex-Gaussian',
    'Shifted inverse Gaussian',
    'Wald race',
    'Two-boundary wiener',
    'Normal',
    'Gamma',
    'Weibull',
    'Log-normal',
    'Beta'
  )
  available_dist[,2] = c(
    'emg',
    'sig',
    'wr',
    'wp',
    'n',
    'ga',
    'we',
    'ln',
    'b'
  )

  # Ex-Gaussian
  common_names = c( 'emg', 'EMG', 'demg', 'pemg',
                    'qemg', 'remg', 'ex-Gaussian',
                    'ex-gaussian', 'ex guass',
                    'ex-gauss', 'exGaussian',
                    'exgaussian' )
  if ( dist %in% common_names ) {
    out = 'emg'
  }

  # Shifted inverse Gaussian
  common_names = c( 'sig', 'SIG', 'dinvgauss', 'pinvgauss',
                    'qinvgauss', 'rinvgauss', 'wald',
                    'Wald' )
  if ( dist %in% common_names ) {
    out = 'sig'
  }

  # Wald race
  common_names = c( 'wr', 'WR', 'dwaldrace', 'pwaldrace',
                    'qwaldrace', 'rwaldrace', 'wald race',
                    'Wald race', 'waldrace' )
  if ( dist %in% common_names ) {
    out = 'wr'
  }

  # Two boundary wiener
  common_names = c( 'wp', 'WP', 'dwiener', 'pwiener',
                    'qwiener', 'rwiener', 'wiener',
                    'Wiener' )
  if ( dist %in% common_names ) {
    out = 'wp'
  }

  # Normal
  common_names = c( 'n', 'N', 'dnorm', 'pnorm',
                    'qnorm', 'rnorm', 'Normal',
                    'normal', 'norm', 'Norm' )
  if ( dist %in% common_names ) {
    out = 'n'
  }

  # Gamma
  common_names = c( 'ga', 'Ga', 'GA', 'dgamma', 'pgamma',
                    'qgamma', 'rgamma', 'Gamma',
                    'gamma' )
  if ( dist %in% common_names ) {
    out = 'ga'
  }

  # Weibull
  common_names = c( 'we', 'We', 'WE', 'dweibull', 'pweibull',
                    'qweibull', 'rweibull', 'Weibull',
                    'weibull', 'weib', 'Weib', 'wei', 'Wei' )
  if ( dist %in% common_names ) {
    out = 'we'
  }

  # Log-normal
  common_names = c( 'ln', 'Ln', 'LN', 'dlnorm', 'plnorm',
                    'qlnorm', 'rlnorm', 'Log-normal',
                    'log-normal', 'log-norm', 'Log-norm' )
  if ( dist %in% common_names ) {
    out = 'ln'
  }

  # Beta
  common_names = c( 'n', 'N', 'dnorm', 'pnorm',
                    'qnorm', 'rnorm', 'Normal',
                    'normal', 'norm', 'Norm' )
  if ( dist %in% common_names ) {
    out = 'n'
  }

  if ( is.null( out ) ) {

    err = paste(
      'Distribution name unknown.', '\n',
      'See below for a list of available ', '\n',
      'distributions and their identifying ', '\n',
      'labels:', '\n',
      paste( paste( available_dist[,1], ' - ',
             available_dist[,2], sep = '' ),
             collapse = '\n' ),
      sep = '' )
    stop( err, call. = FALSE )

  }

  return( out )
}

# Lookup - 06
lower_upper = function( int, dat ) {
  # Purpose:
  # Computes the lower and upper bounds for a set
  # of data within a desired increment.
  # Arguments:
  # int - The increment to round down/up to for
  #       the limits
  # dat - A vector of data
  # Returns:
  # A vector with a lower and upper limit.

  ll = int*floor(min(dat)/int)
  ll = c(ll,int*ceiling(max(dat)/int))

  return( ll )
}
