clear; close all; clc;

%% ========== fixed observations ==========
dt = 1/800;   % original time step

% replace with your actual 15 observation frequencies (Hz)
freq_obs = [100 120 140 160 180 200 220 240 260 280 300 320 340 360 380];
omega_obs = 2*pi*freq_obs(:);   % M x 1
M = length(freq_obs);

%% ========== vary tauR only ==========
tauR_list = [1.38, 1.0, 0.8, 0.5, 0.2, 0.1, 0.05,0.01,0.005,0.001];

%%
K_max = floor(tauR_list ./ dt) +1;
nCase = length(tauR_list);

%%
for ic = 1:nCase
    tauR = tauR_list(ic);
    K = K_max(ic);
    k = 0:K-1;

    % A: M x K
    A = exp(-1i * (omega_obs * (k*dt)));

    s = svd(A, 'econ');
    s_rel = s / s(1);
    cumE = cumsum(s.^2) / sum(s.^2);

    figure;
    plot(cumE,"-o");
    xlabel('Mode index i');
    ylabel('Cum sum of energy');
    title('SVD of A with fixed 15 observations, varying Kmax via \tau_R');
    grid on;
%     legend(label_store, 'Location', 'best');
end



%%


% rel_tol = 1e-3;
% energy_tol = 0.99;
% 
% nCase = length(tauR_list);
% K_all = zeros(nCase,1);
% cond_all = zeros(nCase,1);
% r_dom_all = zeros(nCase,1);
% r_energy_all = zeros(nCase,1);
% rank_all = zeros(nCase,1);
% sigma_store = cell(nCase,1);
% label_store = cell(nCase,1);
% 

% 
%     K_all(ic) = K;
%     cond_all(ic) = s(1)/s(end);
%     r_dom_all(ic) = sum(s_rel > rel_tol);
%     r_energy_all(ic) = find(cumE >= energy_tol, 1, 'first');
%     rank_all(ic) = rank(A);
% 
%     sigma_store{ic} = s_rel;
%     label_store{ic} = sprintf('\\tau_R = %.3f s, K = %d', tauR, K);
% end
% 
% %% ========== table ==========
% T = table(tauR_list(:), K_all, rank_all, r_dom_all, r_energy_all, cond_all, ...
%     'VariableNames', {'tauR_s','K','rankA','r_dom_relTol','r_99energy','condA'});
% disp(T)
% 
% %% ========== plot singular values ==========
% figure;
% hold on;
% for ic = 1:nCase
%     semilogy(1:length(sigma_store{ic}), sigma_store{ic}, '-o', 'LineWidth', 1.5);
% end
% xlabel('Mode index i');
% ylabel('\sigma_i / \sigma_1');
% title('SVD of A with fixed 15 observations, varying K via \tau_R');
% grid on;
% legend(label_store, 'Location', 'best');
% 
% %% ========== plot dominant mode count ==========
% figure;
% plot(K_all, r_dom_all, '-o', 'LineWidth', 1.5);
% xlabel('K');
% ylabel(sprintf('Dominant modes (\\sigma_i/\\sigma_1 > %.0e)', rel_tol));
% title('Dominant SVD modes vs K (fixed 15 observations)');
% grid on;
% 
% %% ========== plot condition number ==========
% figure;
% semilogy(K_all, cond_all, '-o', 'LineWidth', 1.5);
% xlabel('K');
% ylabel('cond(A)');
% title('Condition number vs K (fixed 15 observations)');
% grid on;