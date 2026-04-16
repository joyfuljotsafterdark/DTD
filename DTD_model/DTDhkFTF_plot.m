function DTDhkFTF_plot(out, markerOmegaTauR)
% DTDHKFTF_PLOT
% Part 2: plotting only
%
% Input:
%   out : struct returned by DTDhkFTF_forward
%   markerOmegaTauR : optional marker location on omega*tauR axis

    if nargin < 2
        markerOmegaTauR = 7;
    end

    [~, markerIdx] = min(abs(out.omegaTauR - markerOmegaTauR));

    figure('Color','w','Position',[100 100 1200 500]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    % Left: impulse response
    nexttile([3 1]);
    stem(out.k * out.dt / out.tauR, out.h, 'k', 'filled', ...
        'LineWidth', 1.2, 'MarkerSize', 5);
    xlabel('\tau / \tau_R');
    ylabel('h_k');
    title('Impulse response');
    xlim([0, max(out.k * out.dt / out.tauR)]);

    hmin = min(real(out.h));
    hmax = max(real(out.h));
    if abs(hmax - hmin) < 1e-12
        ylim([hmin-0.5, hmax+0.5]);
    else
        pad = 0.15 * (hmax - hmin);
        ylim([hmin-pad, hmax+pad]);
    end
    grid on;
    box on;

    % Center: Nyquist
%     nexttile([2 1]);
%     plot(out.ReF, out.ImF, 'k', 'LineWidth', 1.5);
%     hold on;
%     plot(out.ReF(markerIdx), out.ImF(markerIdx), 'ko', 'MarkerFaceColor', 'k');
%     hold off;
%     xlabel('Re(F)');
%     ylabel('Im(F)');
%     title('Nyquist plot');
%     axis equal;
%     grid on;
%     box on;

    % Right top: gain
    nexttile;
    plot(out.omegaTauR, out.gain, 'k', 'LineWidth', 1.5);
    hold on;
%     plot(out.omegaTauR(markerIdx), out.gain(markerIdx), 'bo', 'MarkerFaceColor', 'k');
%     hold off;
    ylabel('|F|');
    title('Bode plot');
    xlim([0 50]);
    ylim([-0.5, max(out.gain)+0.5]);
    grid on;
    box on;

    % Right center: wrap phase
    nexttile;
    plot(out.omegaTauR, out.phaseWrapped, 'k', 'LineWidth', 1.5);
    hold on;
%     plot(out.omegaTauR(markerIdx), out.phaseWrapped(markerIdx), 'ko', 'MarkerFaceColor', 'k');
%     hold off;
    ylabel('arg(F)');
    xlim([0 50]);
    ylim([-pi pi]);
    yticks([-pi 0 pi]);
    yticklabels({'-\pi','0','\pi'});
    grid on;
    box on;

     % Right bottom: unwrap phase
    nexttile;
    plot(out.omegaTauR, out.phaseUnwrapped, 'k', 'LineWidth', 1.5);
    hold on;
%     plot(out.omegaTauR(markerIdx), out.phaseWrapped(markerIdx), 'ko', 'MarkerFaceColor', 'k');
%     hold off;
    xlabel('\omega \tau_R');
    ylabel('arg(F)');
    xlim([0 50]);
%     ylim([-pi pi]);
%     yticks([-pi 0 pi]);
%     yticklabels({'-\pi','0','\pi'});
    grid on;
    box on;
end