function result = DTD_hkFTFplot3(amps, idx, k_max)
%DTD_hkFTFplot3
% Three parameters:
%   amps: amplitudes of nonzero h_k terms (sign included)
%   idx: locations k of nonzero h_k terms
%   k_max: maximum index of h_k, so total length is k_max+1
%
% Fixed assumptions in this this version
%   tauR = 1
%   dt   = 0.1 * tauR
%
% Example:
%   amps = [0.4 0.6 0.8 0.6 0.4 -0.1 -0.3 -0.5 -0.3 -0.1];
%   idx  = [1   2   3   4   5    6    7    8    9    10];
%   out = DTD_hkFTFplot3(amps, idx, 10);

    %% -------------------- fixed settings --------------------%%
    tauR = 1;
    dt   = 0.1 * tauR;

    omegaTauR = linspace(0, 25, 3000);
    omega     = omegaTauR ./ tauR;

    % marker point for highlighting (optional, fixed here)
    markerOmegaTauR = 7;   % you can change later
    % --------------------------------------------------------

    %% input check
    if nargin < 3
        error('Usage: result = DTD_hkFTFplot3(amps, idx, k_max)');
    end

    if ~isvector(amps) || ~isvector(idx)
        error('amps and idx must be vectors.');
    end

    amps = amps(:).';   % force row vector
    idx  = idx(:).';    % force row vector

    if length(amps) ~= length(idx)
        error('amps and idx must have the same length.');
    end

    if any(idx < 0) || any(idx > k_max) || any(mod(idx,1) ~= 0)
        error('All idx must be integers in the range [0, k_max].');
    end

    %% build h_k
    k = 0:k_max;
    h = zeros(size(k));

    % accumulate in case repeated indices are provided
    for n = 1:length(idx)
        h(idx(n)+1) = h(idx(n)+1) + amps(n);
    end

    %% compute FTF from full h_k
    % F(omega) = sum_k h_k exp(-i * omega * k * dt)
    E  = exp(-1i * (k(:) * dt) * omega);   % size = (k_max+1) x Nomega
    FTF = h * E;                            % row vector

    gain = abs(FTF);
    phaseWrapped  = angle(FTF);
%     phaseUnwrapped = unwrap(phaseWrapped);

    ReF = real(FTF);
    ImF = imag(FTF);

    %% choose marker index
    [~, markerIdx] = min(abs(omegaTauR - markerOmegaTauR));

    %% plotting
    figure('Color','w','Position',[100 100 1200 500]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % ----- Left: impulse response -----
    nexttile([2 1]);
    stem(k*dt/tauR, h, 'k', 'filled', 'LineWidth', 1.2, 'MarkerSize', 5);
    xlabel('\tau / \tau_R');
    ylabel('h_k');
    title('Impulse response');
    xlim([0, max(k*dt/tauR)]);
%     ylim([min(h)-0.5, max(h)+0.5])

    hmin = min(h);
    hmax = max(h);
    if abs(hmax - hmin) < 1e-12
        ylim([hmin-0.5, hmax+0.5]);
    else
        pad = 0.15 * (hmax - hmin);
        ylim([hmin-pad, hmax+pad]);
    end
    grid on;
    box on;

    % ----- Center: Nyquist -----
%     nexttile([2 1]);
%     plot(ReF, ImF, 'k', 'LineWidth', 1.5);
%     hold on;
%     plot(ReF(markerIdx), ImF(markerIdx), 'bo', 'MarkerFaceColor', 'k');
%     hold off;
%     xlabel('Re(F)');
%     ylabel('Im(F)');
%     title('Nyquist plot');
%     axis equal;
%     grid on;
%     box on;

    % ----- Right top: gain -----
    nexttile;
    plot(omegaTauR, gain, 'k', 'LineWidth', 1.5);
    hold on;
    plot(omegaTauR(markerIdx), gain(markerIdx), 'bo', 'MarkerFaceColor', 'k');
    hold off;
    ylabel('|F|');
    title('Bode plot');
    
    xlim([0 25]);
    ylim([-0.5, max(gain)+0.5])
%     yticks([0 1 2]);
%     yticklabels({'0','1','2'});
    grid on;
    box on;

    % ----- Right bottom: unwrapped phase -----
    nexttile;
    plot(omegaTauR, phaseWrapped, 'k', 'LineWidth', 1.5);
    hold on;
    plot(omegaTauR(markerIdx), phaseWrapped(markerIdx), 'bo', 'MarkerFaceColor', 'k');
    hold off;
    xlabel('\omega \tau_R');
    ylabel('arg(F)');
    xlim([0 25]);
    ylim([-pi pi]);
    yticks([-pi 0 pi]);
    yticklabels({'-\pi','0','\pi'});
    grid on;
    box on;

    %% return results
    result = struct();
    result.tauR = tauR;
    result.dt   = dt;
    result.k    = k;
    result.h    = h;
    result.omegaTauR = omegaTauR;
    result.omega = omega;
    result.FTF   = FTF;
    result.gain  = gain;
    result.phaseWrapped   = phaseWrapped;
    % result.phaseUnwrapped = phaseUnwrapped;
    result.ReF = ReF;
    result.ImF = ImF;
    result.markerIdx = markerIdx;
    result.markerOmegaTauR = omegaTauR(markerIdx);
end