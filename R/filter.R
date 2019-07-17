filter = function(b, a, x, w){
  # Parameters
  # ----------
  #   b : matrix [J, 1]
  #
  #   a : matrix [K, 1]

  #   x : matrix [L, 1]
  #     raw data
  #   w : matrix [M, 1]

  # Returns
  # -------
  #   y : matrix [L, 1]
  #      filtered data
  # Notes
  # -----
  # ORIGINAL DOCUMENTATION FROM THE filter.m SCRIPT
  # Filter a vector.
  # y = filter(b,a,x) returns the solution to the following linear,
  # time-invariant difference equation:
  #
  #    N                   M
  #   sum a(k+1) y(n-k) + sum b(k+1) x(n-k) = 0  for 1<=n<=length(x)
  #   k=0                 k=0
  #
  # where N=length(a)-1 and M=length(b)-1. An equivalent form of this
  # equation is:
  #
  #           N                   M
  #   y(n) = sum c(k+1) y(n-k) + sum d(k+1) x(n-k)  for 1<=n<=length(x)
  #          k=1                 k=0
  #
  # where c = a/a(1) and d = b/a(1).
  #
  # In terms of the z-transform, y is the result of passing the discrete-
  # time signal x through a system characterized by the following rational
  # system function:
  #
  #              M
  #             sum d(k+1) z^(-k)
  #             k=0
  #   H(z) = ----------------------
  #                N
  #           1 + sum c(k+1) z(-k)
  #               k=1
  #
  # [y, sf] = filter(b,a,x,si) sets the initial state of the system, si,
  # and returns the final state, sf.  The state vector is a column vector
  # whose length is equal to the length of the longest coefficient vector
  # minus one.  If si is not set, the initial state vector is set to all
  # zeros.
  #
  # The particular algorithm employed is known as a transposed Direct
  # Form II implementation.
  #
  # SEE ALSO: poly, roots, conv, deconv, residue, polyval, polyderiv, polyinteg

  # Written by Tony Richardson <amr@mpl.ucsd.edu> June 1994.

  # Bug fix by FL (Friedrich.Leisch@ci.tuwien.ac.at) on Oct 12, 1994
  # Ported to R by Dan Broman July 16, 2019
  # ADD ERROR CHECKING
  #   NUMBER OF ARGUMENTS
  #   TYPE OF VAR IN a, b, AND x

  N = length(a)
  M = length(b)
  L = length(x)

  MN = max(c(N, M))
  lw = MN - 1

  # REPLACES postpad. ASSUME postpad is only used to pad and not to trim vector
  # e.g. MN >= length(b)
  if(length(b) < MN){
    b = c(b, rep(0, MN - length(b)))
  }

  # ADD ERROR CHECKING FOR INPUTS? a AND b NEED TO BE COLUMN VECTORS OR Nx1 MATRICIES

  # ADD CHECKING FOR DIMENSIONS OF w? w NEEDS TO BE A COLUMN VECTOR OR Nx1 MATRIX
  # If w is not provided, set initial states to zero
  if(missing(w)){
    # Set initial state to zero
    w = matrix(0, lw, 1)
  }

  # Allocate space for result.
  y = matrix(0, L, 1)

  norm = a[1]
  if(norm == 0){
    stop("First element in second argument must be non-zero.")
  }

  if(norm != 1){
    b = b / norm
  }

  if(N > 1){
    # IIR filter.
    if(length(a) < MN){
      a = c(a, rep(0, MN - length(a)))
    }
    if(norm != 1){
      a = a/norm
    }
    for(index in 1:L){
      y[index, 1] = w[1,] + b[1,] * x[index,]
      # Update state vector
      if(lw > 1){
        w[1:(lw-1), ] = w[2:lw, ] - a[2:lw, ] * y[index, ] + b[2:lw, ] * x[index, ]
        w[lw, ] = b[MN, ] * x[index, ] - a[MN, ] * y[index, ]
      } else {
        w[1, ] = b[MN, ] * x[index, ] - a[MN, ] * y[index, ]
      }
    }
  } else {
    # FIR filter
    if(l1 > 0){
      for(index in 1:L){
        y[index, ] = w[1, ] + b[1, ] * x[index, ]
        w[lw, ] = b[MN] * x[index, ]
      }
    } else {
      y = b[1] * x
    }
  }


  return(y) # ORIGINAL FUNCTION RETURNED y AND w; OCTAVE BASE CODE ONLY RETURNS y
}
