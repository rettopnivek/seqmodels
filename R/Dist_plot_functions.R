#---------------------------------#
# Distribution plotting functions #
#---------------------------------#

# Index
# Lookup - 01:  dist_emg
# Lookup - 02:  dist_sig
# Lookup - 03:  dist_wr
# Lookup - 04:  dist_wp
# Lookup - 05:  dist_n
# Lookup - 06:  dist_ga
# Lookup - 07:  dist_we
# Lookup - 08:  dist_ln
# Lookup - 09:  dist_b
# Lookup - 10:  dist_exp

# Lookup - 01
dist_emg = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the exponentially modified gaussian.
  # Arguments:
  # prm  - A vector of 0-3 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( m = .439, s = .077, l = 2.916 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'emg'
  out$ch = 1

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = demg(
                       x = x,
                       mu = p[1], sigma = p[2], lambda = p[3] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pemg(
        q = x,
        mu = p[1], sigma = p[2], lambda = p[3] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qemg( p = x,
                mu = p[1], sigma = p[2], lambda = p[3] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = demg( x = x, mu = p[1], sigma = p[2], lambda = p[3] )
    p = pemg( q = x, mu = p[1], sigma = p[2], lambda = p[3] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 02
dist_sig = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the shifted inverse gaussian (Wald) distribution.
  # Arguments:
  # prm  - A vector of 0-4 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( k = 0.8, x = 2, t = 0.2, s = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'sig'
  out$ch = 1

  # Truncate x to be within support
  x = x[ x >= 0 ]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dinvgauss(
                       t = x,
                       kappa = p[1], xi = p[2],
                       tau = p[3], sigma = p[4] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list( x = x,
                     y = pinvgauss(
                       t = x,
                       kappa = p[1], xi = p[2],
                       tau = p[3], sigma = p[4] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list( y = x,
                     x = qinvgauss(
                       p = x,
                       kappa = p[1], xi = p[2],
                       tau = p[3], sigma = p[4] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dinvgauss( t = x, kappa = p[1], xi = p[2],
                   tau = p[3], sigma = p[4] )
    p = pinvgauss( t = x, kappa = p[1], xi = p[2],
                   tau = p[3], sigma = p[4] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 03
# FORTHCOMING

# Lookup - 04
dist_wp = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the two boundary wiener process.
  # Arguments:
  # prm  - A vector of 0-5 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the joint density,
  # distribution, quantile, and conditional hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( a = 1.2, z = 0.5, v = 1.0, t0 = 0.3, s = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'wp'
  out$ch = 2
  out$p1 = pwiener( Inf,
                    ch = 1,
                    alpha = p[1],
                    theta = p[2],
                    xi = p[3],
                    tau = p[4],
                    sigma = p[5] )

  # Truncate x to be within support
  x = x[ x >= 0 ]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y1 = dwiener( x, ch = 1,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] ),
                     y0 = dwiener( x, ch = 0,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list( x = x,
                     y1 = pwiener( x, ch = 1,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] ),
                     y0 = pwiener( x, ch = 0,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list( y = x,
                     x1 = qwiener( x, ch = 1,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] ),
                     x0 = qwiener( x, ch = 0,
                                   alpha = p[1],
                                   theta = p[2],
                                   xi = p[3],
                                   tau = p[4],
                                   sigma = p[5] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    k = pwiener( Inf, ch = 1,
                 alpha = p[1],
                 theta = p[2],
                 xi = p[3],
                 tau = p[4],
                 sigma = p[5] )
    d1 = dwiener( x, ch = 1,
                  alpha = p[1],
                  theta = p[2],
                  xi = p[3],
                  tau = p[4],
                  sigma = p[5] )/k
    d0 = dwiener( x, ch = 0,
                  alpha = p[1],
                  theta = p[2],
                  xi = p[3],
                  tau = p[4],
                  sigma = p[5] )/(1-k)
    p1 = pwiener( x, ch = 1,
                  alpha = p[1],
                  theta = p[2],
                  xi = p[3],
                  tau = p[4],
                  sigma = p[5] )/k
    p0 = pwiener( x, ch = 0,
                  alpha = p[1],
                  theta = p[2],
                  xi = p[3],
                  tau = p[4],
                  sigma = p[5] )/(1-k)
    out$h_pv = list( x = x,
                     y1 = d1/(1-p1),
                     y0 = d0/(1-p0)
    )
  }

  return( out )
}

# Lookup - 05
dist_n = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the normal distribution.
  # Arguments:
  # prm  - A vector of 0-2 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( m = 0.0, s = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'n'
  out$ch = 1

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dnorm(
                       x = x,
                       m = p[1], sd = p[2] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pnorm(
        q = x,
        m = p[1], sd = p[2] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qnorm(
        p = x,
        m = p[1], sd = p[2] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dnorm(
      x = x,
      m = p[1], sd = p[2] )
    p = pnorm(
      q = x,
      m = p[1], sd = p[2] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 06
dist_ga = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the gamma distribution.
  # Arguments:
  # prm  - A vector of 0-2 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( s = 1.0, r = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'n'
  out$ch = 1

  # Truncate x to be within support
  x = x[ x >= 0 ]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dgamma(
                       x = x,
                       shape = p[1], rate = p[2] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pgamma(
        q = x,
        shape = p[1], rate = p[2] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qgamma(
        p = x,
        shape = p[1], rate = p[2] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dgamma(
      x = x,
      shape = p[1], rate = p[2] )
    p = pgamma(
      q = x,
      shape = p[1], rate = p[2] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 07
dist_we = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the weibull distribution.
  # Arguments:
  # prm  - A vector of 0-2 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( shape = 1.0, scale = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'n'
  out$ch = 1

  # Truncate x to be within support
  x = x[ x >= 0 ]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dweibull(
                       x = x,
                       shape = p[1], scale = p[2] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pweibull(
        q = x,
        shape = p[1], scale = p[2] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qweibull(
        p = x,
        shape = p[1], scale = p[2] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dweibull(
      x = x,
      shape = p[1], scale = p[2] )
    p = pweibull(
      q = x,
      shape = p[1], scale = p[2] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 08
dist_ln = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the log-normal distribution.
  # Arguments:
  # prm  - A vector of 0-2 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( meanlog = 0.0, sdlog = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'n'
  out$ch = 1

  # Truncate x to be within support
  x = x[ x >= 0 ]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dlnorm(
                       x = x,
                       meanlog = p[1], sdlog = p[2] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = plnorm(
        q = x,
        meanlog = p[1], sdlog = p[2] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qlnorm(
        p = x,
        meanlog = p[1], sdlog = p[2] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dlnorm(
      x = x,
      meanlog = p[1], sdlog = p[2] )
    p = plnorm(
      q = x,
      meanlog = p[1], sdlog = p[2] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 09
dist_b = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the beta distribution.
  # Arguments:
  # prm  - A vector of 0-2 parameters
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( shape1 = 1.0, shape2 = 1.0 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'b'
  out$ch = 1

  # Truncate x to be within support
  x = x[ x >= 0 & x <= 1]

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dbeta(
                       x = x,
                       shape1 = p[1], shape2 = p[2] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pbeta(
        q = x,
        shape1 = p[1], shape2 = p[2] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qbeta(
        p = x,
        shape1 = p[1], shape2 = p[2] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dbeta(
      x = x,
      shape1 = p[1], shape2 = p[2] )
    p = pbeta(
      q = x,
      shape1 = p[1], shape2 = p[2] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

# Lookup - 10
dist_exp = function( prm, type, x, b ) {
  # Purpose:
  # Computes a specified function over a range of values
  # for the exponential function.
  # Arguments:
  # prm  - A vector of 1 parameter
  # type - The type of function to compute
  # x    - An optional vector of values over which to compute
  #        the function
  # b    - An optional vector giving the range of values over
  #        which to compute the function
  # Returns:
  # A list with the x and y-axis values for the density,
  # distribution, quantile, and hazard functions.

  # Initialize output
  out = list(
    d_pv = NULL,
    p_pv = NULL,
    q_pv = NULL,
    h_pv = NULL
  )

  # Define default values
  dv = c( r = 1 )

  p = prm_create( dv, prm )
  ver = check_version( type )
  x = create_x( x, b, ver )
  out$type = ver
  out$prm = p
  out$dist = 'exp'
  out$ch = 1

  # Density
  if ( ver == 1 ) {
    out$d_pv = list( x = x,
                     y = dexp(
                       x = x,
                       rate = p[1] )
    )
  }
  # Distribution function
  if ( ver == 2 ) {
    out$p_pv = list(
      x = x,
      y = pexp(
        q = x,
        rate = p[1] )
    )
  }
  # Quantile function
  if ( ver == 3 ) {
    out$q_pv = list(
      y = x,
      x = qexp(
        p = x,
        rate = p[1] )
    )
  }
  # Hazard function
  if ( ver == 4 ) {
    d = dexp( x = x, rate = p[1] )
    p = pexp( q = x, rate = p[1] )
    out$h_pv = list(
      x = x,
      y = d/(1-p) )
  }

  return( out )
}

