## Copyright (C) 1996, 1997 John W. Eaton
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.

## -*- texinfo -*-
## @deftypefn {Function File} {} fftfilt (@var{b}, @var{x}, @var{n})
##
## With two arguments, @code{fftfilt} filters @var{x} with the FIR filter
## @var{b} using the FFT.
##
## Given the optional third argument, @var{n}, @code{fftfilt} uses the
## overlap-add method to filter @var{x} with @var{b} using an N-point FFT.
##
## If @var{x} is a matrix, filter each column of the matrix.
## @} # deftypefn

## Author: Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
## Created: 3 September 1994
## Adapted-By: jwe


FftFilter <- function(b, n) {
  res = list(b = b, n = n)
  class(res) = "FftFilter"
  res
}

filter.FftFilter <- function(filt, x, ...)
  fftfilt(filt$b, x, filt$n)

# CHANGE
#   x -> a
#   n -> x

fftfilt  <- function(b, a, x = NULL) {

  ## If n is not specified explicitly, we do not use the overlap-add
  ## method at all because loops are really slow.  Otherwise, we only
  ## ensure that the number of points in the FFT is the smallest power
  ## of two larger than N and length(b).  This could result in length
  ## one blocks, but if the user knows better ...
  X = x

  lx = length(x)
  lb = length(b)
  la = length(a)
  n = max(lb, la)
  lrefl = 3 * (n - 1)

  kdc = sum(b) / sum(a)

  if (abs(kdc) < Inf){
    si = rev(cumsum(rev(b - kdc * a)))
  } else {
    si = matrix(0, nrow(a), ncol(a))
  }
  si = si[2:length(si)]

  for (c = 1:ncol(x)){ # filter all columns, one by one
    v = c(2 * x[1, c] - x[seq((lrefl + 1), 2, by = -1), c], x[ ,c],
      2 * x[nrow(x),c] - x[seq((nrow(x) - 1), (nrow(x) - lrefl), by = -1), c]) # a column vector
    ## Do forward and reverse filtering
    v = filter(b, a, v, si * v[1])
  }
## -------------------------------------------2019-07-16
    ## Do forward and reverse filtering
    v = filter(b,a,v,si*v(1));                   # forward filter
    v = flipud(filter(b,a,flipud(v),si*v(end))); # reverse filter
    y(:,c) = v((lrefl+1):(lx+lrefl));




  if (is.null(x)) {
    ## Use FFT with the smallest power of 2 which is >= length (x) +
    ## length (b) - 1 as number of points ...
    N = 2 ^ (ceiling(log(la + lb - 1) / log(2)))
    B = fft(postpad(b, N))
    y = ifft(fft(postpad(a, X)) * B)
  } else {
    ## Use overlap-add method ...
    if (length(x) > 1)
      stop("fftfilt: x has to be a scalar")

    X = 2 ^ (ceiling(log(max(X, lb)) / log(2)))
    L = X - lb + 1
    B = fft(postpad(b, X))
#    B = B[,rep(1., l_x)]
    R = ceiling(la / L)
    y = numeric(la)
    for (r  in  1:R) {
      lo = (r - 1) * L + 1
      hi = min(r * L, la)
      tmp = numeric(0)
      tmp[1:(hi-lo+1)] = a[lo:hi]
      tmp = ifft(fft(postpad(tmp, X)) * B)
      hi  = min(lo+N-1, la)
      y[lo:hi] = y[lo:hi] + tmp[1:(hi-lo+1)]
    }
  }
  y = y[1:la]

  ## Final cleanups: if both a and b are real respectively integer, y
  ## should also be

  if (is.numeric(b) && is.numeric(a))
    y = Re(y)
  if (!any(as.logical(b - round(b)))) {
    ida = !any(as.logical(a - round(a)))
    y[ida] = round(y[ida])
  }
  y
}
