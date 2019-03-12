
function pdf = fb_gampdf (x, a, b)

  if (nargin ~= 3)
    print_usage ;
  end

  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('gampdf: x, a and b must be of common size or scalars');
    end
  end

  sz = size(x);
  pdf = zeros (sz);

  k = find (~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    pdf (k) = NaN;
  end

  k = find ((x > 0) & (a > 0) & (a <= 1) & (b > 0));
  if (any (k))
    if (isscalar(a) && isscalar(b))
      pdf(k) = (x(k) .^ (a - 1)) ...
                .* exp(- x(k) ./ b) ./ gamma (a) ./ (b .^ a);
    else
      pdf(k) = (x(k) .^ (a(k) - 1)) ...
                .* exp(- x(k) ./ b(k)) ./ gamma (a(k)) ./ (b(k) .^ a(k));
    end
  end

  k = find ((x > 0) & (a > 1) & (b > 0));
  if (any (k))
    if (isscalar(a) && isscalar(b))
      pdf(k) = exp (- a .* log (b) + (a-1) .* log (x(k)) ...
                    - x(k) ./ b - gammaln (a));
    else
      pdf(k) = exp (- a(k) .* log (b(k)) + (a(k)-1) .* log (x(k)) ...
                    - x(k) ./ b(k) - gammaln (a(k)));
    end
  end


