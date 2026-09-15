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
colors = [
      0 114 178
    230 159   0
    213  94   0
    204 121 167
      0 158 115
    107  76 154
]/255;
darkText = [17 24 39]/255;
gridTone = [209 213 219]/255;

figureOne = figure( ...
    'Color', 'w', ...
    'InvertHardcopy', 'off', ...
    'Position', [100 100 1120 620]);
axesOne = axes( ...
    'Parent', figureOne, ...
    'Color', 'w', ...
    'XColor', darkText, ...
    'YColor', darkText, ...
    'GridColor', gridTone, ...
    'LineWidth', 0.9, ...
    'FontSize', 10);
hold(axesOne, 'on')

patch(axesOne, ...
    [0 4 4 0], ...
    [settlingBand settlingBand -settlingBand -settlingBand], ...
    [229 231 235]/255, ...
    'FaceAlpha', 0.65, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

for i = 1:numberOfCases
    plot(axesOne, t, xHistory(:, i)/x0, ...
        'LineWidth', 1.7, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s: \\zeta = %.3f', ...
            strrep(caseNames(i), '_', ' '), zetas(i)));
end

yline(axesOne, settlingBand, '--', ...
    'Color', [75 85 99]/255, 'HandleVisibility', 'off');
yline(axesOne, -settlingBand, '--', ...
    'Color', [75 85 99]/255, 'HandleVisibility', 'off');
yline(axesOne, 0, ':', ...
    'Color', [107 114 128]/255, 'HandleVisibility', 'off');

grid(axesOne, 'on')
box(axesOne, 'on')
xlim(axesOne, [0 4])
ylim(axesOne, [-0.82 1.05])
xlabel(axesOne, 'Time, t (s)', 'Color', darkText)
ylabel(axesOne, 'Normalized displacement, x/x_0', 'Color', darkText)
title(axesOne, 'Free response across damping regimes', ...
    'Color', darkText, 'FontWeight', 'bold')
legendOne = legend(axesOne, 'Location', 'eastoutside');
legendOne.Color = 'w';
legendOne.TextColor = darkText;
legendOne.EdgeColor = [156 163 175]/255;
exportgraphics(figureOne, ...
    fullfile(figuresDirectory, 'free_response_comparison.png'), ...
    'Resolution', 300, ...
    'BackgroundColor', 'white')

%% Figure 2: normalized mechanical energy
initialEnergy = 0.5*m*v0^2 + 0.5*k*x0^2;
figureTwo = figure( ...
    'Color', 'w', ...
    'InvertHardcopy', 'off', ...
    'Position', [120 120 1120 620]);
axesTwo = axes( ...
    'Parent', figureTwo, ...
    'Color', 'w', ...
    'XColor', darkText, ...
    'YColor', darkText, ...
    'GridColor', gridTone, ...
    'LineWidth', 0.9, ...
    'FontSize', 10, ...
    'YScale', 'log');
hold(axesTwo, 'on')

for i = 1:numberOfCases
    normalizedEnergy = max(energyHistory(:, i)/initialEnergy, realmin);
    semilogy(axesTwo, t, normalizedEnergy, ...
        'LineWidth', 1.7, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s: \\zeta = %.3f', ...
            strrep(caseNames(i), '_', ' '), zetas(i)));
end

axesTwo.YScale = 'log';
grid(axesTwo, 'on')
box(axesTwo, 'on')
xlim(axesTwo, [0 4])
ylim(axesTwo, [1e-12 1.25])
xlabel(axesTwo, 'Time, t (s)', 'Color', darkText)
ylabel(axesTwo, 'Normalized mechanical energy, E/E_0', 'Color', darkText)
title(axesTwo, 'Mechanical-energy dissipation', ...
    'Color', darkText, 'FontWeight', 'bold')
legendTwo = legend(axesTwo, 'Location', 'eastoutside');
legendTwo.Color = 'w';
legendTwo.TextColor = darkText;
legendTwo.EdgeColor = [156 163 175]/255;
exportgraphics(figureTwo, ...
    fullfile(figuresDirectory, 'energy_decay.png'), ...
    'Resolution', 300, ...
    'BackgroundColor', 'white')

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

figureThree = figure( ...
    'Color', 'w', ...
    'InvertHardcopy', 'off', ...
    'Position', [140 140 1020 620]);
axesThree = axes( ...
    'Parent', figureThree, ...
    'Color', 'w', ...
    'XColor', darkText, ...
    'YColor', darkText, ...
    'GridColor', gridTone, ...
    'LineWidth', 0.9, ...
    'FontSize', 10);
hold(axesThree, 'on')
sweepLine = plot(axesThree, zetaSweep, settlingSweep, ...
    'Color', [0.08 0.35 0.65], ...
    'LineWidth', 1.8);
designMarker = plot(axesThree, zetaDesign, designSettlingTime, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.85 0.24 0.16], ...
    'MarkerEdgeColor', darkText);
xline(axesThree, 1, '--', ...
    'Color', [55 65 81]/255, ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
text(axesThree, 1.025, 13.4, 'Critical damping, \zeta = 1', ...
    'Color', darkText, ...
    'VerticalAlignment', 'top');

grid(axesThree, 'on')
box(axesThree, 'on')
xlim(axesThree, [zetaSweep(1) zetaSweep(end)])
ylim(axesThree, [0 14])
xlabel(axesThree, 'Damping ratio, \zeta', 'Color', darkText)
ylabel(axesThree, 'Two-percent settling time (s)', 'Color', darkText)
title(axesThree, ...
    'Settling-time tradeoff for x(0) = x_0 and v(0) = 0', ...
    'Color', darkText, 'FontWeight', 'bold')
legendThree = legend(axesThree, ...
    [sweepLine designMarker], ...
    {'Numerical sweep', '2% reverse-peak boundary'}, ...
    'Location', 'northeast');
legendThree.Color = 'w';
legendThree.TextColor = darkText;
legendThree.EdgeColor = [156 163 175]/255;
exportgraphics(figureThree, ...
    fullfile(figuresDirectory, 'settling_time_sweep.png'), ...
    'Resolution', 300, ...
    'BackgroundColor', 'white')

savefig(figureOne, fullfile(figuresDirectory, 'free_response_comparison.fig'));

savefig(figureTwo, fullfile(figuresDirectory, 'energy_decay.fig'));

savefig(figureThree, fullfile(figuresDirectory, 'settling_time_sweep.fig'));

% Capture the numerical environment and verification errors for this run.
verification = table(caseNames, ode45MaxError, ...
    'VariableNames', {'case_name', 'max_displacement_error_m'});
writetable(verification, fullfile(dataDirectory, 'ode45_verification.csv'));
fid = fopen(fullfile(dataDirectory, 'run_environment.txt'), 'w');
assert(fid ~= -1, 'Could not open environment output.');
fprintf(fid, 'MATLAB %s\nRelease %s\nPlatform %s\n', version, version('-release'), computer);
fprintf(fid, 'ode45 RelTol=1e-10, AbsTol=1e-12; tVerify=linspace(0,4,2001)\n');
fclose(fid);
