function ts = compute_settling_time(t, x, tolerance, referenceMagnitude)
%COMPUTE_SETTLING_TIME Earliest permanent entry into a tolerance band.
%
%   ts = COMPUTE_SETTLING_TIME(t,x,tolerance,referenceMagnitude) returns the
%   earliest time after which abs(x) remains less than or equal to
%
%       tolerance * abs(referenceMagnitude).
%
%   The function returns NaN when the response has not settled before the
%   final supplied time. A linear interpolation estimates the final band
%   crossing between adjacent samples.


validateattributes(t, {'numeric'}, {'vector','real','finite','nonnegative'}, ...
    mfilename, 't');
validateattributes(x, {'numeric'}, {'vector','real','finite'}, ...
    mfilename, 'x');
validateattributes(tolerance, {'numeric'}, ...
    {'scalar','real','finite','positive','<',1}, mfilename, 'tolerance');
validateattributes(referenceMagnitude, {'numeric'}, ...
    {'scalar','real','finite','nonzero'}, mfilename, 'referenceMagnitude');


if numel(t) ~= numel(x)
    error('compute_settling_time:SizeMismatch', ...
        't and x must contain the same number of samples.');
end


t = t(:);
x = x(:);


if any(diff(t) <= 0)
    error('compute_settling_time:NonMonotonicTime', ...
        't must be strictly increasing.');
end


band = tolerance*abs(referenceMagnitude);
distanceFromBand = abs(x) - band;
lastOutside = find(distanceFromBand > 0, 1, 'last');


if isempty(lastOutside)
    ts = t(1);
    return
end


if lastOutside == numel(t)
    ts = NaN;
    return
end


t1 = t(lastOutside);
t2 = t(lastOutside + 1);
y1 = distanceFromBand(lastOutside);
y2 = distanceFromBand(lastOutside + 1);


% y1 is positive and y2 is non-positive at the final band crossing.
if y1 == y2
    ts = t2;
else
    ts = t1 - y1*(t2 - t1)/(y2 - y1);
end
end