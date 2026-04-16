clear all;
close all;
clc;

%% ========== load data ============ %%

raw_data = readmatrix("data\test3_FTF_output.csv");

freq = raw_data(3:end,2);

%% ========== load data ============ %%
omega = 2*pi*freq;
% dt = 1./(2.*freq);
dt = 5e-5;
%tauR = 1.38;
tauR = 0.015;
k_max = floor(tauR/dt) +1;

%%
for n=1:15
    A = zeros(n,k_max(n));
    k = 0:k_max(n)-1;
    for m=1:k_max(n)
        A(:,m) = exp(omega(1:n).*(-i)*dt(n)*m);
    end
    [u,sig,vh] = svd(A,"econ");
    S = diag(sig).^2;
    temp_sum = cumsum(S)./sum(S);
    figure();
    plot(temp_sum);
end

%%

for m=1:k_max
    A(:,m) = exp(omega(:).*(-i)*dt*m);
end

%%
[u,sig,vh] = svd(A,"econ");

S = diag(sig);

figure();
plot(S)

%%
temp_sum = cumsum(S)./sum(S);

plot(temp_sum);



