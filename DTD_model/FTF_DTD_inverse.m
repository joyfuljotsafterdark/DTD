clear;
close all;
clc;

%% =========================================================
% 1. LOAD DATA
% ==========================================================
ftfMatrix = readmatrix("./data/test3_FTF_output.csv");

% Data columns
frequencyHz       = ftfMatrix(3:end, 2);
ftfGainMeasured   = ftfMatrix(3:end, 3);
ftfPhaseDeg       = ftfMatrix(3:end, 4);

% Optional: wrap degree phase into [-180, 180]
% ftfPhaseDeg = mod(ftfPhaseDeg + 180, 360) - 180;

% Convert degree to rad
ftfPhaseWrapped   = deg2rad(ftfPhaseDeg);
ftfPhaseUnwrapped = unwrap(ftfPhaseWrapped);

% Construct measured complex FTF
ftfMeasured = ftfGainMeasured .* exp(1i * ftfPhaseWrapped);

%% =========================================================
% 2. SETTINGS
% ==========================================================
useSumConstraint        = false;
timeStep                = 5e-5;     % s
regularizationLambda    = 1e-4;
numberOfLowFreqPoints   = 3;
dtdLength               = 10;
sumConstraintWeight     = 1e2;

angularFrequency        = 2 * pi * frequencyHz;
numberOfFrequencyPoints = length(frequencyHz);

%% =========================================================
% 3. ESTIMATE BULK DELAY
% phase ≈ intercept - omega * bulkDelay
% ==========================================================
lowFreqOmega = angularFrequency(1:numberOfLowFreqPoints);
lowFreqPhase = ftfPhaseUnwrapped(1:numberOfLowFreqPoints);

phaseFitCoeff = polyfit(lowFreqOmega, lowFreqPhase, 1);
phaseSlope    = phaseFitCoeff(1);
phaseIntercept = phaseFitCoeff(2); 

bulkDelay = -phaseSlope;

fprintf('Estimated bulk delay = %.6e s\n', bulkDelay);

bulkDelayResponse = exp(-1i * angularFrequency * bulkDelay);

%% =========================================================
% 4. REMOVE BULK DELAY
% residualMeasured = F / exp(-i*omega*tau0)
% ==========================================================
residualMeasured      = ftfMeasured ./ bulkDelayResponse;
residualGainMeasured  = abs(residualMeasured);
residualPhaseMeasured = unwrap(angle(residualMeasured));

%% =========================================================
% 5. BUILD LINEAR SYSTEM
% residualMeasured(omega_j) = sum_{k=0}^{K-1} residualDtd(k) * exp(-i*omega_j*k*dt)
% ==========================================================
designMatrix = zeros(numberOfFrequencyPoints, dtdLength);

for freqIndex = 1:numberOfFrequencyPoints
    for delayIndex = 0:dtdLength-1
        designMatrix(freqIndex, delayIndex+1) = ...
            exp(-1i * angularFrequency(freqIndex) * delayIndex * timeStep);
    end
end

% Second-difference smoothing matrix
if dtdLength >= 3
    smoothMatrix = zeros(dtdLength-2, dtdLength);
    for rowIndex = 1:dtdLength-2
        smoothMatrix(rowIndex, rowIndex:rowIndex+2) = [1, -2, 1];
    end
else
    smoothMatrix = zeros(0, dtdLength);
end

augmentedMatrix = [
    designMatrix;
    sqrt(regularizationLambda) * smoothMatrix
];

augmentedRhs = [
    residualMeasured;
    zeros(size(smoothMatrix, 1), 1)
];

if useSumConstraint
    augmentedMatrix = [
        augmentedMatrix;
        sumConstraintWeight * ones(1, dtdLength)
    ];

    augmentedRhs = [
        augmentedRhs;
        sumConstraintWeight * 1
    ];
end

% Convert complex least squares to real least squares
realSystemMatrix = [
    real(augmentedMatrix);
    imag(augmentedMatrix)
];

realSystemRhs = [
    real(augmentedRhs);
    imag(augmentedRhs)
];

%% =========================================================
% 6. SOLVE FOR RESIDUAL DTD
% ==========================================================
residualDtd = realSystemMatrix \ realSystemRhs;

%% =========================================================
% 7. RECONSTRUCT RESPONSE
% ==========================================================
residualFitted = designMatrix * residualDtd;
ftfFitted      = bulkDelayResponse .* residualFitted;

ftfGainFitted      = abs(ftfFitted);
ftfPhaseFitted     = unwrap(angle(ftfFitted));
bulkDelayPhaseOnly = unwrap(angle(bulkDelayResponse));

%% =========================================================
% 8. ERROR REPORT
% ==========================================================
relativeErrorFtf      = norm(ftfFitted - ftfMeasured) / norm(ftfMeasured);
relativeErrorResidual = norm(residualFitted - residualMeasured) / norm(residualMeasured);

fprintf('Relative error in FTF      = %.6e\n', relativeErrorFtf);
fprintf('Relative error in residual = %.6e\n', relativeErrorResidual);
fprintf('Sum of residual DTD        = %.6f\n', sum(residualDtd));

disp('Estimated residual DTD coefficients = ');
disp(residualDtd(:).');

%% =========================================================
% 9. PLOTS
% ==========================================================
figure('Color', 'w', 'Position', [80 80 1300 900]);

% -------- Original FTF: Gain --------
subplot(3,2,1);
plot(frequencyHz, ftfGainMeasured, 'ko', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(frequencyHz, ftfGainFitted, 'b-', 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('|F|');
title('Original FTF: Gain');
legend('Measured', 'Fitted', 'Location', 'best');
grid on;

% -------- Original FTF: Phase --------
subplot(3,2,2);
plot(frequencyHz, ftfPhaseUnwrapped, 'ko', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(frequencyHz, ftfPhaseFitted, 'b-', 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('Phase [rad]');
title('Original FTF: Phase');
legend('Measured', 'Fitted', 'Location', 'best');
grid on;

% -------- Bulk delay only --------
subplot(3,2,3);
plot(frequencyHz, ftfPhaseUnwrapped, 'ko', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(frequencyHz, bulkDelayPhaseOnly, 'r--', 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('Phase [rad]');
title(sprintf('Bulk Delay Only, tau_0 = %.4e s', bulkDelay));
legend('Measured phase', 'Pure delay phase', 'Location', 'best');
grid on;

% -------- Residual response --------
subplot(3,2,4);
yyaxis left;
plot(frequencyHz, residualGainMeasured, 'ko', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(frequencyHz, abs(residualFitted), 'b-', 'LineWidth', 1.8);
ylabel('|G|');

yyaxis right;
plot(frequencyHz, residualPhaseMeasured, 'rs', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(frequencyHz, unwrap(angle(residualFitted)), 'm-', 'LineWidth', 1.6);
ylabel('Phase(G) [rad]');

xlabel('f [Hz]');
title('Residual Response after Removing Bulk Delay');
legend('|G| measured', '|G| fitted', 'phase(G) measured', 'phase(G) fitted', 'Location', 'best');
grid on;

% -------- Residual DTD --------
subplot(3,2,5);
stem(0:dtdLength-1, residualDtd, 'filled', 'LineWidth', 1.2);
xlabel('k');
ylabel('g_k');
title(sprintf('Estimated Residual DTD, K = %d', dtdLength));
grid on;

% -------- Nyquist --------
subplot(3,2,6);
plot(real(ftfMeasured), imag(ftfMeasured), 'ko', 'LineWidth', 1.2, 'MarkerSize', 7); hold on;
plot(real(ftfFitted), imag(ftfFitted), 'b-', 'LineWidth', 1.8);
xlabel('Re(F)');
ylabel('Im(F)');
title('Nyquist: Original FTF');
legend('Measured', 'Fitted', 'Location', 'best');
grid on;
axis equal;

%% =========================================================
% 10. OPTIONAL: RESIDUAL NYQUIST
% ==========================================================
figure('Color', 'w', 'Position', [150 150 600 500]);
plot(real(residualMeasured), imag(residualMeasured), 'ko', 'LineWidth', 1.2, 'MarkerSize', 7); hold on;
plot(real(residualFitted), imag(residualFitted), 'b-', 'LineWidth', 1.8);
xlabel('Re(G)');
ylabel('Im(G)');
title('Nyquist of Residual Response');
legend('Measured residual', 'Fitted residual', 'Location', 'best');
grid on;
axis equal;

%% =========================================================
% 11. OPTIONAL: PHYSICAL TIME AXIS
% ==========================================================
residualDelayTime = (0:dtdLength-1).' * timeStep;

figure('Color', 'w', 'Position', [220 220 650 450]);
stem(residualDelayTime * 1000, residualDtd, 'filled', 'LineWidth', 1.2);
xlabel('Residual delay time [ms]');
ylabel('g_k');
title('Residual DTD on Physical Time Axis');
grid on;

