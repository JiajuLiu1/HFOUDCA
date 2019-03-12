
function cdf = fb_gamcdf (x, a, b)

  if (nargin ~= 3)
    print_usage ;
  end

  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('gamcdf: x, a and b must be of common size or scalars');
    end
  end

  sz = size (x);
  cdf = zeros (sz);

  k = find (~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    cdf (k) = NaN;
  end

  k = find ((x > 0) & (a > 0) & (b > 0));
  if (any (k))
    if (isscalar (a) && isscalar(b))
      cdf (k) = gammainc (x(k) ./ b, a);
    else
      cdf (k) = gammainc (x(k) ./ b(k), a(k));
    end
  end

