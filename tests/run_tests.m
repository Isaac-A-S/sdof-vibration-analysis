%% Validation checks for the analytical SDOF response
% Run from any folder with:
%   run("tests/run_tests.m")




clear
clc




testPath = mfilename('fullpath');
testsDirectory = fileparts(testPath);
projectRoot = fileparts(testsDirectory);
addpath(fullfile(projectRoot, 'matlab'))




m = 1000;
k = 40000;
x0 = 0.50;
v0 = 0;
wn = sqrt(k/m);
cCritical = 2*sqrt(k*m);




cCases = [0, 2000, 6000, 12000, cCritical, 16000];
t = linspace(0, 4, 4001).';
odeOptions = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);




fprintf('Running SDOF response checks...\n')




%% Check 1: initial conditions and governing acceleration
for c = cCases
    [x, v, a] = free_response(t, m, k, c, x0, v0);




    assert(abs(x(1) - x0) < 1e-12, ...
        sprintf('Initial displacement check failed for c = %.6g.', c))
    assert(abs(v(1) - v0) < 1e-11, ...
        sprintf('Initial velocity check failed for c = %.6g.', c))




    expectedInitialAcceleration = -(c*v0 + k*x0)/m;
    assert(abs(a(1) - expectedInitialAcceleration) < 1e-10, ...
        sprintf('Initial acceleration check failed for c = %.6g.', c))
end




fprintf('  PASS: initial conditions and initial acceleration\n')




%% Check 2: analytical response versus ODE45
for c = cCases
    stateDerivative = @(~, state) [
        state(2)
        -(c/m)*state(2) - (k/m)*state(1)
    ];




    [~, numericalState] = ode45( ...
        stateDerivative, t, [x0; v0], odeOptions);




    xExact = free_response(t, m, k, c, x0, v0);
    maximumError = max(abs(xExact - numericalState(:, 1)));




    assert(maximumError < 1e-8, ...
        sprintf(['Analytical/ODE45 mismatch for c = %.6g. ' ...
        'Maximum error was %.3e m.'], c, maximumError))
end




fprintf('  PASS: analytical response agrees with ODE45\n')




%% Check 3: mechanical energy never increases
initialEnergy = 0.5*m*v0^2 + 0.5*k*x0^2;
energyTolerance = 1e-10*initialEnergy;




for c = cCases
    [x, v] = free_response(t, m, k, c, x0, v0);
    energy = 0.5*m.*v.^2 + 0.5*k.*x.^2;




    assert(max(diff(energy)) <= energyTolerance, ...
        sprintf('Mechanical energy increased for c = %.6g.', c))
end




fprintf('  PASS: mechanical energy is non-increasing\n')




%% Check 4: solution continuity near critical damping
criticalOffset = 1e-6;
tCritical = linspace(0, 2, 2001).';




xBelow = free_response( ...
    tCritical, m, k, cCritical*(1 - criticalOffset), x0, v0);
xAt = free_response( ...
    tCritical, m, k, cCritical, x0, v0);
xAbove = free_response( ...
    tCritical, m, k, cCritical*(1 + criticalOffset), x0, v0);




relativeBelowDifference = max(abs(xBelow - xAt))/abs(x0);
relativeAboveDifference = max(abs(xAbove - xAt))/abs(x0);




assert(relativeBelowDifference < 5e-6, ...
    'Underdamped branch is not continuous near critical damping.')
assert(relativeAboveDifference < 5e-6, ...
    'Overdamped branch is not continuous near critical damping.')




fprintf('  PASS: response branches are continuous near zeta = 1\n')




%% Check 5: two-percent damping design
allowedReversePeak = 0.02;
logPeak = -log(allowedReversePeak);
zetaDesign = logPeak/sqrt(pi^2 + logPeak^2);
cDesign = zetaDesign*cCritical;




tFine = (0:1e-4:4).';
xDesign = free_response(tFine, m, k, cDesign, x0, v0);
settlingTime = compute_settling_time(tFine, xDesign, 0.02, x0);
reversePeakRatio = exp(-pi*zetaDesign/sqrt(1 - zetaDesign^2));




assert(abs(reversePeakRatio - allowedReversePeak) < 1e-12, ...
    'The analytical reverse-peak design is incorrect.')
assert(abs(settlingTime - 0.5696) < 5e-4, ...
    'The two-percent settling-time regression check failed.')




fprintf('  PASS: two-percent reverse-peak design\n')




%% Check 6: overdamped response is monotonic for the selected conditions
xOverdamped = free_response(t, m, k, 16000, x0, v0);




assert(all(xOverdamped >= -1e-12), ...
    'Selected overdamped response crossed equilibrium.')
assert(all(diff(xOverdamped) <= 1e-12), ...
    'Selected overdamped response was not monotonic.')




fprintf('  PASS: selected overdamped response is monotonic\n')
fprintf('All SDOF response checks passed.\n')
