function [x, v, a, meta] = free_response(t, m, k, c, x0, v0)
%FREE_RESPONSE Linear SDOF free response: m*x'' + c*x' + k*x = 0.
% t [s] may be a scalar or array; outputs retain its shape.
% m > 0 [kg], k > 0 [N/m], c >= 0 [N*s/m]. Initial conditions:
% x(0) = x0 [m], v(0) = v0 [m/s]. No add-on toolbox is required.
% Inputs are evaluated in double precision. Extreme parameter magnitudes
% that overflow derived quantities are outside the intended domain.

validateattributes(t, {'numeric'}, {'real','finite','nonnegative','nonempty'}, mfilename, 't');
validateattributes(m, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'm');
validateattributes(k, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'k');
validateattributes(c, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'c');
validateattributes(x0, {'numeric'}, {'scalar','real','finite'}, mfilename, 'x0');
validateattributes(v0, {'numeric'}, {'scalar','real','finite'}, mfilename, 'v0');
t = double(t); m = double(m); k = double(k); c = double(c);
x0 = double(x0); v0 = double(v0);
wn = sqrt(k/m);
cCritical = 2*sqrt(k*m);
zeta = c/cCritical;
alpha = c/(2*m);

% Only a few floating-point spacings use the repeated-root limit.
criticalTolerance = 8*eps(1);
if zeta < 1 - criticalTolerance
    wd = wn*sqrt((1-zeta)*(1+zeta));
    u = wd*t;
    sincu = ones(size(u));
    small = abs(u) < 1e-4;
    sincu(small) = 1-u(small).^2/6+u(small).^4/120;
    sincu(~small) = sin(u(~small))./u(~small);
    S = t.*sincu; % sin(wd*t)/wd, with its continuous wd -> 0 limit
    decay = exp(-alpha*t);
    x = decay.*(x0*cos(u)+(v0+alpha*x0)*S);
    v = decay.*(v0*cos(u)-(alpha*v0+wn^2*x0)*S);
    regime = "underdamped";
elseif zeta > 1 + criticalTolerance
    q = sqrt((zeta-1)*(zeta+1));
    rSlow = -wn/(zeta+q); % rationalized root avoids cancellation at large zeta
    rFast = -wn*(zeta+q);
    gap = 2*wn*q;
    slow = exp(rSlow*t);
    D = slow.*(-expm1(-gap*t))/gap; % difference of exponentials / root gap
    B = v0-rSlow*x0;
    x = x0*slow+B*D;
    v = rSlow*x+B*exp(rFast*t);
    wd = NaN;
    regime = "overdamped";
else
    B = v0+wn*x0;
    decay = exp(-wn*t);
    x = (x0+B*t).*decay;
    v = (v0-wn*B*t).*decay;
    wd = NaN;
    regime = "critical";
end

% Acceleration comes from the ODE, not independent differentiation.
% Checking this same identity alone would therefore be a circular test.
a = -(c*v+k*x)/m;
meta = struct('wn',wn,'fn',wn/(2*pi),'cCritical',cCritical, ...
    'zeta',zeta,'wd',wd,'regime',regime);
end
