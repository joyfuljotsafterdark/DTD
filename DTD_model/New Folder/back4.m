clear; close all; clc;

%% ================= user inputs =================
tauR = 1.38;                  % keep tauR fixed
f_band = 140;                 % fixed comparison band (Hz), safer than 150
rel_tol = 1e-3;               % dominant mode threshold: sigma_i/sigma_1 > rel_tol
energy_tol = 0.99;            % cumulative energy threshold

% Your observation frequencies in Hz
% Example:
freq_obs = [100 120 140 160 180 200 220 240 260 280 300 320 340 360 380];
% Replace with your actual frequencies

% Candidate effective dt values
dt_list = [5e-5, 1/800, 1/500, 1/300];

%% ================= fixed low-frequency observations =================
freq_use = freq_obs(freq_obs <= f_band);
omega_use = 2*pi*freq_use(:);     % column vector
M = length(freq_use);

if M == 0
    error('No observation frequencies are inside the selected low-frequency band.');
end

fprintf('Using %d observations in %.1f Hz band:\n', M, f_band);
disp(freq_use)

%% ================= loop over dt_eff =================
nCase = length(dt_list);

K_all        = zeros(nCase,1);
cond_all     = zeros(nCase,1);
r_dom_all    = zeros(nCase,1);   % dominant mode count by relative threshold
r_energy_all = zeros(nCase,1);   % mode count reaching 99%% energy
rank_all     = zeros(nCase,1);   % numerical rank from MATLAB

sigma_store = cell(nCase,1);
label_store = cell(nCase,1);

for ic = 1:nCase
    dt_eff = dt_list(ic);

    % Nyquist check for the fixed band
    fN = 1/(2*dt_eff);
    if fN < f_band
        warning('dt_eff = %.6g s has Nyquist %.2f Hz < f_band %.2f Hz. Skip.', ...
                dt_eff, fN, f_band);
        K_all(ic)        = NaN;
        cond_all(ic)     = NaN;
        r_dom_all(ic)    = NaN;
        r_energy_all(ic) = NaN;
        rank_all(ic)     = NaN;
        sigma_store{ic}  = [];
        label_store{ic}  = sprintf('\\Deltat=%.4g s (invalid)', dt_eff);
        continue;
    end

    % K = floor(tauR/dt_eff)+1
    K_eff = floor(tauR/dt_eff) + 1;
    k = 0:K_eff-1;

    % Build A: M x K
    A = exp(-1i * (omega_use * (k*dt_eff)));

    % SVD
    s = svd(A, 'econ');
    s_rel = s / s(1);
    cumE = cumsum(s.^2) / sum(s.^2);

    % Metrics
    K_all(ic) = K_eff;
    cond_all(ic) = s(1) / s(end);
    r_dom_all(ic) = sum(s_rel > rel_tol);
    r_energy_all(ic) = find(cumE >= energy_tol, 1, 'first');
    rank_all(ic) = rank(A);

    sigma_store{ic} = s_rel;
    label_store{ic} = sprintf('\\Deltat=%.4g s, K=%d', dt_eff, K_eff);
end

%% ================= results table =================
T = table(dt_list(:), K_all, rank_all, r_dom_all, r_energy_all, cond_all, ...
    'VariableNames', {'dt_eff_s','K_eff','rankA','r_dom_relTol','r_99energy','condA'});
disp(T)

%% ================= plot 1: singular value decay =================
figure;
hold on;
for ic = 1:nCase
    if ~isempty(sigma_store{ic})
        semilogy(1:length(sigma_store{ic}), sigma_store{ic}, '-o', 'LineWidth', 1.5);
    end
end
xlabel('Mode index i');
ylabel('\sigma_i / \sigma_1');
title(sprintf('Normalized singular values of A (fixed band: f <= %.1f Hz)', f_band));
grid on;
legend(label_store, 'Location', 'best');

%% ================= plot 2: dominant mode count vs K =================
figure;
plot(K_all, r_dom_all, '-o', 'LineWidth', 1.5);
xlabel('K_{eff}');
ylabel(sprintf('Dominant modes (\\sigma_i/\\sigma_1 > %.0e)', rel_tol));
title('Dominant mode count vs K_{eff}');
grid on;

%% ================= plot 3: condition number vs K =================
figure;
semilogy(K_all, cond_all, '-o', 'LineWidth', 1.5);
xlabel('K_{eff}');
ylabel('cond(A)');
title('Condition number vs K_{eff}');
grid on;

%% ================= optional: cumulative energy plot =================
figure;
hold on;
for ic = 1:nCase
    if ~isempty(sigma_store{ic})
        s_rel = sigma_store{ic};
        % Recover absolute shape not needed; use normalized values squared
        cumE_rel = cumsum(s_rel.^2) / sum(s_rel.^2);
        plot(1:length(cumE_rel), cumE_rel, '-o', 'LineWidth', 1.5);
    end
end
xlabel('Mode index i');
ylabel('Cumulative energy');
title('Cumulative singular value energy');
grid on;
legend(label_store, 'Location', 'best');