filtfilt = function(b, a, x) {
  # Parameters
  # ----------
  #   b : matrix [J, 1]
  #
  #   a : matrix [K, 1]

  #   x : matrix [L, 1]
  #     raw data

  # Returns
  # -------
  #   y : matrix [L, 1]
  #      filtered data
  # Notes
  # -----
  # ORIGINAL DOCUMENTATION FROM THE filtfilt.m SCRIPT
  # Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
  # Copyright (C) 2007 Francesco Potortì <pot@gnu.org>
  # Copyright (C) 2008 Luca Citi <lciti@essex.ac.uk>
  #
  # This program is free software: you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.
  #
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  # GNU General Public License for more details.
  #
  # You should have received a copy of the GNU General Public License
  # along with this program; see the file COPYING. If not, see
  # <https://www.gnu.org/licenses/>.

  # -*- texinfo -*-
  # @deftypefn {Function File} {@var{y} =} filtfilt (@var{b}, @var{a}, @var{x})
  #
  # Forward and reverse filter the signal. This corrects for phase
  # distortion introduced by a one-pass filter, though it does square the
  # magnitude response in the process. That's the theory at least.  In
  # practice the phase correction is not perfect, and magnitude response
  # is distorted, particularly in the stop band.
  #
  # Example
  # @example
  # @group
  # [b, a]=butter(3, 0.1);                  # 5 Hz low-pass filter
  # t = 0:0.01:1.0;                         # 1 second sample
  # x=sin(2*pi*t*2.3)+0.25*randn(size(t));  # 2.3 Hz sinusoid+noise
  # y = filtfilt(b,a,x); z = filter(b,a,x); # apply filter
  # plot(t,x,';data;',t,y,';filtfilt;',t,z,';filter;')
  # @end group
  # @end example
  # @end deftypefn
  # Ported to R by Dan Broman July 17, 2019

  lx = length(x)
  lb = length(b)
  la = length(a)
  n = max(lb, la)
  lrefl = 3 * (n - 1)
  if(la < n){
    a = c(a, rep(0, n - la))
  }
  if(lb < n){
    b = c(b, rep(0, n - lb))
  }

  # Allocate space for result.
  y = matrix(0, lx, ncol(x))

  ## Compute a the initial state taking inspiration from
  ## Likhterov & Kopeika, 2003. "Hardware-efficient technique for
  ##     minimizing startup transients in Direct Form II digital filters"
  kdc = sum(b) / sum(a)
  if (abs(kdc) < Inf){ # neither NaN nor +/- Inf
    si = rev(cumsum(rev(b - kdc * a)))
  } else {
    si = matrix(0, nrow(a), ncol(a))
  }
  si = si[2:length(si)]

  for(c in 1:ncol(x)){ # filter all columns, one by one
    v = c(2 * x[1, c] - x[seq((lrefl + 1), 2, by = -1), c], x[ ,c],
      2 * x[nrow(x),c] - x[seq((nrow(x) - 1), (nrow(x) - lrefl), by = -1), c]) # a column vector

    ## Do forward and reverse filtering
    v = filter(b, a, matrix(v), matrix(si * v[1]))                    # forward filter
    v = rev(filter(b, a, matrix(rev(v)), matrix(si * v[length(v)])))  # reverse filter
    y[ , c] = v[(lrefl + 1):(lx + lrefl)]
  }

  return(y)
}
