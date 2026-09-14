%% Transient response and damping design of an SDOF oscillator
% This script compares analytical free responses across damping regimes,
% verifies them against ODE45, computes engineering response metrics, and
% exports the summary data and figures used by the repository.


clear
close all
clc


%% Project paths
scriptPath = mfilename('fullpath');
matlabDirectory = fileparts(scriptPath);
projectRoot = fileparts(matlabDirectory);
dataDirectory = fullfile(projectRoot, 'data');
figuresDirectory = fullfile(projectRoot, 'figures');


addpath(matlabDirectory)


if ~isfolder(dataDirectory)
    mkdir(dataDirectory)
end
if ~isfolder(figuresDirectory)
    mkdir(figuresDirectory)
end


%% Illustrative system parameters
m = 1000;          % equivalent mass [kg]
k = 40000;         % equivalent stiffness [N/m]
x0 = 0.50;         % initial displacement [m]
v0 = 0;            % initial velocity [m/s]
settlingBand = 0.02;


wn = sqrt(k/m);
fn = wn/(2*pi);
cCritical = 2*sqrt(k*m);


% If a first reverse peak of 2% is acceptable, this damping ratio places
% that peak exactly on the tolerance boundary.
allowedReversePeak = 0.02;
logPeak = -log(allowedReversePeak);
zetaDesign = logPeak/sqrt(pi^2 + logPeak^2);
cDesign = zetaDesign*cCritical;


caseNames = [
    "low_damping"
    "moderate_damping"
    "two_percent_design"
    "near_critical"
    "critical"
    "overdamped"
];


cValues = [
    2000
    6000
    cDesign
    12000
    cCritical
    16000
];


%% Time bases
% A fine uniform grid makes the settling-time definition reproducible.
t = (0:1e-4:4).';
tVerify = linspace(0, 4, 2001).';


numberOfCases = numel(cValues);
xHistory = zeros(numel(t), numberOfCases);
vHistory = zeros(numel(t), numberOfCases);
energyHistory = zeros(numel(t), numberOfCases);


zetas = zeros(numberOfCases, 1);
regimes = strings(numberOfCases, 1);
dampedFrequencies = NaN(numberOfCases, 1);
envelopeHalfLives = NaN(numberOfCases, 1);
settlingTimes = NaN(numberOfCases, 1);
firstReversePeakPct = NaN(numberOfCases, 1);
ode45MaxError = NaN(numberOfCases, 1);


odeOptions = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);


%% Evaluate each damping case
for i = 1:numberOfCases
    c = cValues(i);


    [x, v, ~, meta] = free_response(t, m, k, c, x0, v0);
    xHistory(:, i) = x;
    vHistory(:, i) = v;


    energy = 0.5*m.*v.^2 + 0.5*k.*x.^2;
    energyHistory(:, i) = energy;


    zetas(i) = meta.zeta;
    regimes(i) = meta.regime;
    dampedFrequencies(i) = meta.wd;
    settlingTimes(i) = compute_settling_time( ...
        t, x, settlingBand, x0);


    if meta.zeta < 1
        envelopeHalfLives(i) = log(2)/(meta.zeta*meta.wn);
        firstReversePeakPct(i) = 100*exp( ...
            -pi*meta.zeta/sqrt(1 - meta.zeta^2));
    end


    % Independent numerical verification. The second-order model is written
    % as y1_dot = y2 and y2_dot = -(c/m)y2 -(k/m)y1.
    stateDerivative = @(~, state) [
        state(2)
        -(c/m)*state(2) - (k/m)*state(1)
    ];


    [~, numericalState] = ode45( ...
        stateDerivative, tVerify, [x0; v0], odeOptions);


    xExact = free_response(tVerify, m, k, c, x0, v0);
    ode45MaxError(i) = max(abs(xExact - numericalState(:, 1)));
end


%% Export case summary
summary = table( ...
    caseNames, ...
    cValues, ...
    zetas, ...
    regimes, ...
    dampedFrequencies, ...
    envelopeHalfLives, ...
    settlingTimes, ...
    firstReversePeakPct, ...
    'VariableNames', { ...
        'case_name', ...
        'c_Ns_per_m', ...
        'zeta', ...
        'regime', ...
        'omega_d_rad_per_s', ...
        'envelope_half_life_s', ...
        'settling_time_2pct_s', ...
        'first_reverse_peak_pct'});


writetable(summary, fullfile(dataDirectory, 'case_summary.csv'))


disp(summary)
fprintf('\nNatural frequency: %.6f rad/s (%.6f Hz)\n', wn, fn)
fprintf('Critical damping: %.6f N*s/m\n', cCritical)
fprintf('Maximum analytical-versus-ODE45 displacement error: %.3e m\n', ...
    max(ode45MaxError))
fprintf(['Two-percent reverse-peak design: zeta = %.6f, ' ...
    'c = %.3f N*s/m, settling time = %.4f s\n'], ...
    zetaDesign, cDesign, settlingTimes(caseNames == "two_percent_design"))


%% Figure 1: normalized displacement
colors = lines(numberOfCases);
figureOne = figure('Color', 'w', 'Position', [100 100 980 580]);
hold on


for i = 1:numberOfCases
    plot(t, xHistory(:, i)/x0, ...
        'LineWidth', 1.7, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s: \\zeta = %.3f', ...
            strrep(caseNames(i), '_', ' '), zetas(i)));
end


yline(settlingBand, '--k', 'HandleVisibility', 'off');
yline(-settlingBand, '--k', 'HandleVisibility', 'off');
yline(0, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');


grid on
box on
xlim([0 4])
xlabel('Time, t (s)')
ylabel('Normalized displacement, x/x_0')
title('Free response across damping regimes')
legend('Location', 'eastoutside')
exportgraphics(figureOne, ...
    fullfile(figuresDirectory, 'free_response_comparison.png'), ...
    'Resolution', 300)


%% Figure 2: normalized mechanical energy
initialEnergy = 0.5*m*v0^2 + 0.5*k*x0^2;
figureTwo = figure('Color', 'w', 'Position', [120 120 980 580]);
hold on


for i = 1:numberOfCases
    normalizedEnergy = max(energyHistory(:, i)/initialEnergy, realmin);
    semilogy(t, normalizedEnergy, ...
        'LineWidth', 1.7, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s: \\zeta = %.3f', ...
            strrep(caseNames(i), '_', ' '), zetas(i)));
end


grid on
box on
xlim([0 4])
xlabel('Time, t (s)')
ylabel('Normalized mechanical energy, E/E_0')
title('Mechanical-energy dissipation')
legend('Location', 'eastoutside')
exportgraphics(figureTwo, ...
    fullfile(figuresDirectory, 'energy_decay.png'), ...
    'Resolution', 300)


%% Figure 3: damping-ratio design sweep
zetaSweep = linspace(0.05, 2.00, 391).';
tSweep = (0:1e-3:15).';
settlingSweep = NaN(size(zetaSweep));


for i = 1:numel(zetaSweep)
    cSweep = zetaSweep(i)*cCritical;
    xSweep = free_response(tSweep, m, k, cSweep, x0, v0);
    settlingSweep(i) = compute_settling_time( ...
        tSweep, xSweep, settlingBand, x0);
end


designResponse = free_response(tSweep, m, k, cDesign, x0, v0);
designSettlingTime = compute_settling_time( ...
    tSweep, designResponse, settlingBand, x0);


figureThree = figure('Color', 'w', 'Position', [140 140 900 560]);
plot(zetaSweep, settlingSweep, ...
    'Color', [0.08 0.35 0.65], ...
    'LineWidth', 1.8)
hold on
plot(zetaDesign, designSettlingTime, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.85 0.24 0.16], ...
    'MarkerEdgeColor', 'k')
xline(1, '--k', 'Critical damping, \zeta = 1', ...
    'LabelVerticalAlignment', 'bottom');


grid on
box on
xlim([zetaSweep(1) zetaSweep(end)])
xlabel('Damping ratio, \zeta')
ylabel('Two-percent settling time (s)')
title('Settling-time tradeoff for x(0) = x_0 and v(0) = 0')
legend('Numerical sweep', '2% reverse-peak boundary', ...
    'Location', 'northeast')
exportgraphics(figureThree, ...
    fullfile(figuresDirectory, 'settling_time_sweep.png'), ...
    'Resolution', 300)