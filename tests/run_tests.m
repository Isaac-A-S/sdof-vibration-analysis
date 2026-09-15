%% Verification checks for the analytical SDOF response
% Run from the repository root with:
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


%% Check 7: general initial states against independent matrix exponential
% expm solves the constant-coefficient state system without using our roots.
A = @(c) [0, 1; -k/m, -c/m];
for zeta = [0, 0.4, 1-1e-10, 1, 1+1e-10, 2, 100]
    for initialState = [0.5, 0, -0.3; -1.2, 2, 0.7]
        for tt = [0, 0.01, 0.3, 2]
            [xx,vv] = free_response(tt,m,k,zeta*cCritical,initialState(1),initialState(2));
            reference = expm(A(zeta*cCritical)*tt)*initialState;
            assert(norm([xx;vv]-reference,inf) < 1e-9, 'State/expm mismatch.');
        end
    end
end
fprintf('  PASS: general initial states and near-critical branches versus expm\n')

%% Check 8: differential consistency and undamped energy conservation
h = 1e-5;
for c = cCases
    tt = linspace(0.1,2,101).';
    [~,vv,aa] = free_response(tt,m,k,c,0.3,-0.8);
    [xp,vp] = free_response(tt+h,m,k,c,0.3,-0.8);
    [xm,vm] = free_response(tt-h,m,k,c,0.3,-0.8);
    assert(max(abs((xp-xm)/(2*h)-vv)) < 1e-7, 'dx/dt does not match v.');
    assert(max(abs((vp-vm)/(2*h)-aa)) < 1e-6, 'dv/dt does not match a.');
end
[xu,vu] = free_response(t,m,k,0,0.3,-0.8);
Eu = 0.5*m*vu.^2+0.5*k*xu.^2;
assert(max(abs(Eu-Eu(1)))/Eu(1) < 1e-12, 'Undamped energy drift.');
assert(isequal(size(free_response([0,1,2],m,k,2000,x0,v0)),[1,3]));
assert(isequal(size(free_response([0;1;2],m,k,2000,x0,v0)),[3,1]));
fprintf('  PASS: differential consistency, conservation, and output shapes\n')

%% Check 9: settling helper behavior, re-entry, and input rejection
assert(compute_settling_time([0,1],[0.01,0],0.02,1) == 0);
assert(isnan(compute_settling_time([0,1],[1,0.1],0.02,1)));
assert(abs(compute_settling_time([0,1,2,3],[1,0,0.1,0],0.02,1)-2.8) < 1e-12);
assert(abs(compute_settling_time([0,1],[-1,0],0.02,1)-0.98) < 1e-12);
rejected = false;
try
    compute_settling_time([0,0],[1,0],0.02,1);
catch exception
    rejected = strcmp(exception.identifier,'compute_settling_time:NonMonotonicTime');
end
assert(rejected, 'Repeated time was not rejected.');
rejected = false;
try
    compute_settling_time([],[],0.02,1);
catch
    rejected = true;
end
assert(rejected, 'Empty record was not rejected.');
fprintf('  PASS: settling helper edge cases\n')

%% Check 10: design peak from solution, grid refinement, and tail bounds
wdDesign = wn*sqrt(1-zetaDesign^2);
tPeak = pi/wdDesign;
[xPeak,vPeak] = free_response(tPeak,m,k,cDesign,x0,0);
assert(abs(xPeak/x0+0.02) < 1e-12 && abs(vPeak) < 1e-10);
% Deliberately include the exact tangential peak in the sample record.
tWithPeak = sort([tFine;tPeak]);
tsPeak = compute_settling_time(tWithPeak,free_response(tWithPeak,m,k,cDesign,x0,0),0.02,x0);
assert(abs(tsPeak-settlingTime) < 1e-6, 'Roundoff changed design settling.');
for c = [cCases(2:end),cDesign]
    previous = NaN;
    for dt = [1e-3,5e-4,1e-4]
        tg = (0:dt:8).';
        xx = free_response(tg,m,k,c,x0,0);
        current = compute_settling_time(tg,xx,0.02,x0);
        if ~isnan(previous)
            assert(abs(current-previous) < 1e-4,'Settling grid did not converge.');
        end
        previous = current;
    end
    short = compute_settling_time(tFine,free_response(tFine,m,k,c,x0,0),0.02,x0);
    assert(abs(short-current) < 1e-6,'Longer horizon changed case settling.');
end
% For every underdamped sweep case this decreasing envelope is already
% below the band at 15 s. Critical/overdamped v0=0 cases are monotone.
for zeta = linspace(0.05,0.995,190)
    assert(exp(-zeta*wn*15)/sqrt(1-zeta^2) < 0.02);
end
fprintf('  PASS: actual design peak, grid refinement, and sampled-case horizons\n')
fprintf('All 10 SDOF response check groups passed.\n')
