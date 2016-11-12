#-------------------#
# Utility functions #
#-------------------#

dist_calc = function(t,ch,prm,dist,prb,ver) {
  # Purpose:
  # Calculates different functions for a distribution given
  # a vector of parameters.
  # Arguments:
  # t    - A sorted sequence of times
  # ch   - A choice value (0 or 1)
  # prm  - A vector of parameters for a given distribution
  # dist - A distribution label
  # prb  - A sequence of quantile probabilities
  # ver  - The type of function to calculate for the distribution
  # Returns:
  # A list of times and the values for the function of the distribution.

  stop_f = T # Check if the supplied label matches a distribution

  # Wald race
  if ( dist == 'wr' ) {

    stop_f = F

    # Define default values
    pv = c( k1 = 1, xi1 = 3, tau1 = 0,
            k0 = 1, xi0 = 1, tau0 = 0,
            s1 = 1, s0 = 1, rl = 0 )
    # Extract inputted parameters
    nms = names( prm )
    if ( length( nms ) == 0 ) {
      nms = names( pv )[1:length(prm)]
      names( prm ) = nms
    }
    pv[ nms ] = prm[ nms ]

    if ( ver == 'CDF' ) {

      y = pwaldrace( t, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )

    }

    if ( ver == 'PDF' ) {

      y = dwaldrace( t, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )

    }

    if ( ver == 'QPE' ) {

      p = pwaldrace( Inf, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )
      t = qwaldrace( prb*p, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )
      y = prb*p
    }

    if ( ver == 'HF' ) {

      p = pwaldrace( Inf, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )
      d = dwaldrace( t, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )
      D = pwaldrace( t, ch,
                     pv[1], pv[2], pv[3], pv[4],
                     pv[5], pv[6], pv[7], pv[8], pv[9] )
      y = (d/p)/(1-D/p)

    }
  }

  # Inverse gaussian
  if ( dist == 'ig' ) {

    stop_f = F

    # Define default values
    pv = c( k = 1, xi = 3, s = 1 )
    # Extract inputted parameters
    nms = names( prm )
    if ( length( nms ) == 0 ) {
      nms = names( pv )[1:length(prm)]
      names( prm ) = nms
    }
    pv[ nms ] = prm[ nms ]

    if ( ver == 'CDF' ) {

      y = pinvgauss( t, pv[1], pv[2], pv[3] )

    }

    if ( ver == 'PDF' ) {

      y = dinvgauss( t, pv[1], pv[2], pv[3] )

    }

    if ( ver == 'QPE' ) {

      t = qinvgauss( prb, pv[1], pv[2], pv[3] )
      y = prb
    }

    if ( ver == 'HF' ) {

      d = dinvgauss( t, pv[1], pv[2], pv[3] )
      D = pinvgauss( t, pv[1], pv[2], pv[3] )
      y = d/(1-D)

    }
  }

  # Exponentially modified gaussian
  if ( dist == 'emg' ) {

    stop_f = F

    # Define default values
    pv = c( mu = -2, s = 1, l = 1 )
    # Extract inputted parameters
    nms = names( prm )
    if ( length( nms ) == 0 ) {
      nms = names( pv )[1:length(prm)]
      names( prm ) = nms
    }
    pv[ nms ] = prm[ nms ]

    if ( ver == 'CDF' ) {

      y = pemg( t, pv[1], pv[2], pv[3] )

    }

    if ( ver == 'PDF' ) {

      y = demg( t, pv[1], pv[2], pv[3] )

    }

    if ( ver == 'QPE' ) {

      t = qemg( prb, pv[1], pv[2], pv[3] )
      y = prb
    }

    if ( ver == 'HF' ) {

      d = demg( t, pv[1], pv[2], pv[3] )
      D = pemg( t, pv[1], pv[2], pv[3] )
      y = d/(1-D)

    }
  }

  if (stop_f) out = list() else out = list( t = t, y = y )
  return( out )
}
