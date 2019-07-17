filtfilt = function(b, a, x) {
  # Parameters
  # ----------
  #   b : filter
  #
  #   a : filter

  #   x : filter
  #     raw data

  # Returns
  # -------
  #   y : filter
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
  y = rep(0, lx)

  ## Compute a the initial state taking inspiration from
  ## Likhterov & Kopeika, 2003. "Hardware-efficient technique for
  ##     minimizing startup transients in Direct Form II digital filters"
  kdc = sum(b) / sum(a)
  if (abs(kdc) < Inf){ # neither NaN nor +/- Inf
    si = rev(cumsum(rev(b - kdc * a)))
  } else {
    si = rep(0, la)
  }
  si = si[2:length(si)]

  # ORIGINAL CODE LOOPED OVER COLUMNS OF x IF DATA WERE LARGER THAN Nx1
  v = c(2 * x[1] - x[seq((lrefl + 1), 2, by = -1)], x,
    2 * x[length(x)] - x[seq((length(x) - 1), (length(x) - lrefl), by = -1)]) # a column vector

  ## Do forward and reverse filtering
  v = filter(b, a, v, si * v[1])                    # forward filter
  v = rev(filter(b, a, rev(v), si * v[length(v)]))  # reverse filter
  y = v[(lrefl + 1):(lx + lrefl)]

  return(y)
}
