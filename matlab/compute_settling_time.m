function ts = compute_settling_time(t, x, tolerance, referenceMagnitude)
%COMPUTE_SETTLING_TIME Last sampled tolerance-band entry, interpolated.
% Band: abs(x) <= tolerance*abs(referenceMagnitude).
% Returns NaN if the last sample is outside; t(1) if all are inside.
% This finite-record estimate cannot exclude excursions between samples or
% after t(end). Check time-grid refinement and the model's tail separately.
% A 32-ULP displacement allowance handles roundoff at a tangential peak;
% it is numerical slack, not physical design margin.

validateattributes(t, {'numeric'}, {'vector','real','finite','nonnegative','nonempty'}, mfilename, 't');
validateattributes(x, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'x');
validateattributes(tolerance, {'numeric'}, {'scalar','real','finite','positive','<',1}, mfilename, 'tolerance');
validateattributes(referenceMagnitude, {'numeric'}, {'scalar','real','finite','nonzero'}, mfilename, 'referenceMagnitude');
if numel(t) ~= numel(x) || numel(t) < 2
    error('compute_settling_time:SizeMismatch', 'Supply at least two matching samples.');
end
 t = double(t(:)); x = double(x(:));
if any(diff(t) <= 0)
    error('compute_settling_time:NonMonotonicTime', 't must be strictly increasing.');
end
referenceMagnitude = abs(double(referenceMagnitude));
band = double(tolerance)*referenceMagnitude;
slack = 32*eps(referenceMagnitude);
lastOutside = find(abs(x) > band+slack, 1, 'last');
if isempty(lastOutside)
    ts = t(1);
elseif lastOutside == numel(t)
    ts = NaN;
else
    j = lastOutside;
    boundary = sign(x(j))*band;
    fraction = (boundary-x(j))/(x(j+1)-x(j));
    fraction = min(1,max(0,fraction));
    ts = t(j)+fraction*(t(j+1)-t(j));
end
end
