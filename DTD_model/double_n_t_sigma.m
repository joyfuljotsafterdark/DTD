clear;
close all;
clc;

%% =========================================================
% 1. LOAD DATA
% ==========================================================
ftfMatrix = readmatrix("./data/test1_FTF_output.csv");

frequencyHz     = ftfMatrix(3:end, 2);
ftfGainMeasured = ftfMatrix(3:end, 3);
ftfPhaseDeg     = ftfMatrix(3:end, 4);

% Optional: wrap phase in degree into [-180, 180]
% ftfPhaseDeg = mod(ftfPhaseDeg + 180, 360) - 180;

ftfPhaseWrapped   = deg2rad(ftfPhaseDeg);
ftfPhaseUnwrapped = unwrap(ftfPhaseWrapped);

ftfMeasured = ftfGainMeasured .* exp(1i * ftfPhaseWrapped);

angularFrequency = 2 * pi * frequencyHz;

%% =========================================================
% 2. SETTINGS
% ==========================================================
useLowFrequencyConstraint = false;   % true / false
lowFrequencyConstraintWeight = 1e1;  % only used if true

% Parameter vector:
% x = [pathAmplitude1, pathDelay1, pathSpread1, ...
%      pathAmplitude2, pathDelay2, pathSpread2]

initialGuess = [ ...
    1.0, 1.0e-3, 3.0e-4, ...
   -0.8, 3.0e-3, 6.0e-4  ...
];

lowerBound = [ ...
   -5.0, 0.0e-3, 1.0e-6, ...
   -5.0, 0.0e-3, 1.0e-6  ...
];

upperBound = [ ...
    5.0, 10.0e-3, 5.0e-3, ...
    5.0, 10.0e-3, 5.0e-3  ...
];

optimizationOptions = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

%% =========================================================
% 3. FIT DOUBLE n-tau-sigma MODEL
% ==========================================================
objectiveFunction = @(parameterVector) ...
    buildDoubleNtaoSigmaResidual( ...
        parameterVector, ...
        angularFrequency, ...
        ftfMeasured, ...
        useLowFrequencyConstraint, ...
        lowFrequencyConstraintWeight);

fittedParameterVector = lsqnonlin( ...
    objectiveFunction, ...
    initialGuess, ...
    lowerBound, ...
    upperBound, ...
    optimizationOptions);

%% =========================================================
% 4. EXTRACT FITTED PARAMETERS
% ==========================================================
pathAmplitude1 = fittedParameterVector(1);
pathDelay1     = fittedParameterVector(2);
pathSpread1    = fittedParameterVector(3);

pathAmplitude2 = fittedParameterVector(4);
pathDelay2     = fittedParameterVector(5);
pathSpread2    = fittedParameterVector(6);

fprintf('\n===== Fitted Parameters =====\n');
fprintf('Path 1 amplitude = %.6f\n', pathAmplitude1);
fprintf('Path 1 delay     = %.6e s\n', pathDelay1);
fprintf('Path 1 spread    = %.6e s\n', pathSpread1);
fprintf('Path 2 amplitude = %.6f\n', pathAmplitude2);
fprintf('Path 2 delay     = %.6e s\n', pathDelay2);
fprintf('Path 2 spread    = %.6e s\n', pathSpread2);

%% =========================================================
% 5. RECONSTRUCT FITTED FTF
% ==========================================================
ftfPath1 = pathAmplitude1 .* exp( ...
    -1i * angularFrequency * pathDelay1 ...
    -0.5 * (angularFrequency.^2) * (pathSpread1^2));

ftfPath2 = pathAmplitude2 .* exp( ...
    -1i * angularFrequency * pathDelay2 ...
    -0.5 * (angularFrequency.^2) * (pathSpread2^2));

ftfFitted = ftfPath1 + ftfPath2;

ftfGainFitted      = abs(ftfFitted);
ftfPhaseFitted     = unwrap(angle(ftfFitted));
ftfPhasePath1      = unwrap(angle(ftfPath1));
ftfPhasePath2      = unwrap(angle(ftfPath2));
ftfGainPath1       = abs(ftfPath1);
ftfGainPath2       = abs(ftfPath2);

relativeErrorFtf = norm(ftfFitted - ftfMeasured) / norm(ftfMeasured);
fprintf('Relative error in FTF = %.6e\n', relativeErrorFtf);

%% =========================================================
% 6. BUILD EQUIVALENT TIME-DOMAIN DISTRIBUTION h(t)
% For F(omega) = a * exp(-i*omega*tau - 0.5*omega^2*sigma^2)
% Corresponding h(t) = a/(sqrt(2*pi)*sigma) * exp(-(t-tau)^2/(2*sigma^2))
% ==========================================================
timeMax = max([pathDelay1 + 4*pathSpread1, pathDelay2 + 4*pathSpread2, 6e-3]);
timeAxis = linspace(0, timeMax, 2000).';

hPath1 = pathAmplitude1 ./ (sqrt(2*pi) * pathSpread1) .* ...
    exp(-(timeAxis - pathDelay1).^2 ./ (2 * pathSpread1^2));

hPath2 = pathAmplitude2 ./ (sqrt(2*pi) * pathSpread2) .* ...
    exp(-(timeAxis - pathDelay2).^2 ./ (2 * pathSpread2^2));

hTotal = hPath1 + hPath2;

%% =========================================================
% 7. PLOTS
% ==========================================================
figure('Color', 'w', 'Position', [80 80 1300 900]);

% -------- Gain --------
subplot(3,2,1);
plot(frequencyHz, ftfGainMeasured, 'ko', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(frequencyHz, ftfGainFitted, 'b-', 'LineWidth', 1.8);
plot(frequencyHz, ftfGainPath1, 'r--', 'LineWidth', 1.2);
plot(frequencyHz, ftfGainPath2, 'm--', 'LineWidth', 1.2);
xlabel('f [Hz]');
ylabel('|F|');
title('Double n-\tau-\sigma Fit: Gain');
legend('Measured', 'Total fit', 'Path 1', 'Path 2', 'Location', 'best');
grid on;

% -------- Phase --------
subplot(3,2,2);
plot(frequencyHz, ftfPhaseUnwrapped, 'ko', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(frequencyHz, ftfPhaseFitted, 'b-', 'LineWidth', 1.8);
plot(frequencyHz, ftfPhasePath1, 'r--', 'LineWidth', 1.2);
plot(frequencyHz, ftfPhasePath2, 'm--', 'LineWidth', 1.2);
xlabel('f [Hz]');
ylabel('Phase [rad]');
title('Double n-\tau-\sigma Fit: Phase');
legend('Measured', 'Total fit', 'Path 1', 'Path 2', 'Location', 'best');
grid on;

% -------- Nyquist --------
subplot(3,2,3);
plot(real(ftfMeasured), imag(ftfMeasured), 'ko', 'LineWidth', 1.2, 'MarkerSize', 7); hold on;
plot(real(ftfFitted), imag(ftfFitted), 'b-', 'LineWidth', 1.8);
plot(real(ftfPath1), imag(ftfPath1), 'r--', 'LineWidth', 1.2);
plot(real(ftfPath2), imag(ftfPath2), 'm--', 'LineWidth', 1.2);
xlabel('Re(F)');
ylabel('Im(F)');
title('Nyquist');
legend('Measured', 'Total fit', 'Path 1', 'Path 2', 'Location', 'best');
grid on;
axis equal;

% -------- Time-domain distribution --------
subplot(3,2,4);
plot(timeAxis * 1000, hTotal, 'b-', 'LineWidth', 1.8); hold on;
plot(timeAxis * 1000, hPath1, 'r--', 'LineWidth', 1.2);
plot(timeAxis * 1000, hPath2, 'm--', 'LineWidth', 1.2);
xlabel('Delay time [ms]');
ylabel('h(t)');
title('Equivalent Time-Delay Distribution');
legend('Total', 'Path 1', 'Path 2', 'Location', 'best');
grid on;

% -------- Real / Imag versus frequency --------
subplot(3,2,5);
plot(frequencyHz, real(ftfMeasured), 'ko-', 'LineWidth', 1.0, 'MarkerSize', 5); hold on;
plot(frequencyHz, real(ftfFitted), 'b-', 'LineWidth', 1.8);
plot(frequencyHz, imag(ftfMeasured), 'rs-', 'LineWidth', 1.0, 'MarkerSize', 5);
plot(frequencyHz, imag(ftfFitted), 'm-', 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('Value');
title('Real / Imag Comparison');
legend('Re measured', 'Re fitted', 'Im measured', 'Im fitted', 'Location', 'best');
grid on;

% -------- Text summary --------
subplot(3,2,6);
axis off;
text(0.00, 0.90, sprintf('Path 1 amplitude = %.4f', pathAmplitude1), 'FontSize', 11);
text(0.00, 0.78, sprintf('Path 1 delay     = %.4e s', pathDelay1), 'FontSize', 11);
text(0.00, 0.66, sprintf('Path 1 spread    = %.4e s', pathSpread1), 'FontSize', 11);

text(0.00, 0.48, sprintf('Path 2 amplitude = %.4f', pathAmplitude2), 'FontSize', 11);
text(0.00, 0.36, sprintf('Path 2 delay     = %.4e s', pathDelay2), 'FontSize', 11);
text(0.00, 0.24, sprintf('Path 2 spread    = %.4e s', pathSpread2), 'FontSize', 11);

text(0.00, 0.06, sprintf('Relative error   = %.4e', relativeErrorFtf), 'FontSize', 11);
title('Fitted Parameters');

%% =========================================================
% LOCAL FUNCTION
% ==========================================================
function residualVector = buildDoubleNtaoSigmaResidual( ...
    parameterVector, ...
    angularFrequency, ...
    ftfMeasured, ...
    useLowFrequencyConstraint, ...
    lowFrequencyConstraintWeight)

    pathAmplitude1 = parameterVector(1);
    pathDelay1     = parameterVector(2);
    pathSpread1    = parameterVector(3);

    pathAmplitude2 = parameterVector(4);
    pathDelay2     = parameterVector(5);
    pathSpread2    = parameterVector(6);

    ftfModel = ...
        pathAmplitude1 .* exp( ...
            -1i * angularFrequency * pathDelay1 ...
            -0.5 * (angularFrequency.^2) * (pathSpread1^2)) ...
        + ...
        pathAmplitude2 .* exp( ...
            -1i * angularFrequency * pathDelay2 ...
            -0.5 * (angularFrequency.^2) * (pathSpread2^2));

    complexResidual = ftfModel - ftfMeasured;

    residualVector = [
        real(complexResidual);
        imag(complexResidual)
    ];

    % Optional low-frequency constraint:
    % F(0) = a1 + a2 ~= 1
    if useLowFrequencyConstraint
        lowFrequencyResidual = lowFrequencyConstraintWeight * ...
            (pathAmplitude1 + pathAmplitude2 - 1);
        residualVector = [
            residualVector;
            lowFrequencyResidual
        ];
    end
end