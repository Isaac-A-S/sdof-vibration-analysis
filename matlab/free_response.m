function [x, v, a, meta] = free_response(t, m, k, c, x0, v0)
%FREE_RESPONSE Closed-form response of a linear mass-spring-damper system.
%
%   [x,v,a,meta] = FREE_RESPONSE(t,m,k,c,x0,v0) evaluates
%
%       m*x_ddot + c*x_dot + k*x = 0
%
%   at the times in t for initial conditions x(0)=x0 and x_dot(0)=v0.
%   The function handles underdamped, critically damped, and overdamped
%   regimes without requiring a MATLAB add-on toolbox.
%
%   Inputs
%       t   - scalar or array of times [s], all t >= 0
%       m   - mass [kg], m > 0
%       k   - stiffness [N/m], k > 0
%       c   - viscous damping coefficient [N*s/m], c >= 0
%       x0  - initial displacement [m]
%       v0  - initial velocity [m/s]
%
%   Outputs
%       x   - displacement [m], same size as t
%       v   - velocity [m/s], same size as t
%       a   - acceleration [m/s^2], same size as t
%       meta - structure containing wn, fn, cCritical, zeta, wd, and regime


validateattributes(t,  {'numeric'}, {'real','finite','nonnegative'}, mfilename, 't');
validateattributes(m,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'm');
validateattributes(k,  {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'k');
validateattributes(c,  {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'c');
validateattributes(x0, {'numeric'}, {'scalar','real','finite'}, mfilename, 'x0');
validateattributes(v0, {'numeric'}, {'scalar','real','finite'}, mfilename, 'v0');


wn = sqrt(k/m);
fn = wn/(2*pi);
cCritical = 2*sqrt(k*m);
zeta = c/cCritical;
alpha = zeta*wn;


% Use a small tolerance so values affected only by roundoff are treated as
% critical. Values outside the tolerance use their exact analytical branch.
criticalTolerance = 1e-8;


if zeta < 1 - criticalTolerance
    wd = wn*sqrt(1 - zeta^2);


    A = x0;
    B = (v0 + alpha*x0)/wd;


    decay = exp(-alpha.*t);
    cosine = cos(wd.*t);
    sine = sin(wd.*t);


    x = decay .* (A.*cosine + B.*sine);
    v = decay .* ( ...
        -alpha.*(A.*cosine + B.*sine) ...
        - A.*wd.*sine + B.*wd.*cosine);


    regime = "underdamped";


elseif zeta > 1 + criticalTolerance
    rootTerm = wn*sqrt(zeta^2 - 1);
    r1 = -alpha + rootTerm;
    r2 = -alpha - rootTerm;


    C1 = (v0 - r2*x0)/(r1 - r2);
    C2 = (r1*x0 - v0)/(r1 - r2);


    x = C1.*exp(r1.*t) + C2.*exp(r2.*t);
    v = C1.*r1.*exp(r1.*t) + C2.*r2.*exp(r2.*t);


    wd = NaN;
    regime = "overdamped";


else
    % Limit of the analytical solution as zeta approaches one.
    A = x0;
    B = v0 + wn*x0;


    decay = exp(-wn.*t);
    x = (A + B.*t).*decay;
    v = (B - wn.*(A + B.*t)).*decay;


    wd = NaN;
    regime = "critical";
end


% Evaluate acceleration from the governing equation. This also avoids a
% second numerical differentiation of the analytical velocity.
a = -(c.*v + k.*x)./m;


meta = struct( ...
    'wn', wn, ...
    'fn', fn, ...
    'cCritical', cCritical, ...
    'zeta', zeta, ...
    'wd', wd, ...
    'regime', regime);
end