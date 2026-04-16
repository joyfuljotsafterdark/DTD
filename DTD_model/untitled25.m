clear all;
clc;
close;

%%
raw_data = readmatrix("data\test3_FTF_output.csv");
FTF_phase = raw_data(3:end,4);
%% ========== fixed observations ==========
dt = 5e-5;   % original time step

% replace with your actual 15 observation frequencies (Hz)
freq_obs = [100 120 140 160 180 200 220 240 260 280 300 320 340 360 380];
omega_obs = 2*pi*freq_obs(:);   % M x 1
M = length(freq_obs);

%%
tauR = 0.01;

K_max = floor(tauR/dt) +1;

niter = 10;
%%
h_recons = zeros(niter,K_max);

for ic = 1:niter
    k = 0:K_max-1;

    % A: M x K
    A = exp(-1i * (omega_obs * (k*dt)));

    h_recon = mldivide(A,FTF_phase);
    h_recon = real(h_recon);

    h_recons(ic,:) = h_recon';
    
end

%%