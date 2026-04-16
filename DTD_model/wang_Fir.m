clear all;
close all;
clc;

%%
k = 1:30;
dt = 0.5;
tau1 = 3.5;
sigma1 = 0.6;
tauc = 7;
taus1 = 1.5;
taus2 = 2.4;

tau2 = taus1 + tauc;
tau3 = taus2 + tauc;

sigma2 = taus1 / 3;
sigma3 = taus2 / 3;

%%

dist1  = (dt / (sigma1*sqrt(2*pi))) .* exp(-(k.*dt-tau1).^2/(2*sigma1^2));
dist2  = (dt / (sigma2*sqrt(2*pi))) .* exp(-(k.*dt-tau2).^2/(2*sigma2^2));
dist3  = -(dt / (sigma3*sqrt(2*pi))) .* exp(-(k.*dt-tau3).^2/(2*sigma3^2));



h = dist1 + dist2 +dist3;

%%
time = k.*dt;

figure();
scatter(time, dist3);

%%

omega = linspace(0,25,1000);

temp_sum = zeros(1,1000);

for k=1:30
    temp_sum = temp_sum + h(k).*exp(-i*k.*omega*dt);
end


%%
gain = abs(temp_sum);
wrap_phase = angle(temp_sum);
unwrap_phase = unwrap(wrap_phase);

%%
figure();
plot(omega,gain);

%%
figure();
plot(omega,wrap_phase);
xlim([0 5]);


%%
figure();
plot(omega,rad2deg(unwrap_phase));
xlim([0 5]);


%%