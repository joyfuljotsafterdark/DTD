clear all;
close all;
clc;
%%
tauR = 1;
dt = 0.1 * tauR;

%%
k_max = 10;
k = 0:k_max;
h = zeros(size(k));
h(k==2) = 1.4;
h(k==6) = -0.4;

%%
omegaTauR = linspace(0,25,2000);
omega = omegaTauR./tauR;
FTF = 1.4 * exp(-1i * omega * (2*dt))+ (-0.4) * exp(-1i * omega * (6*dt));

gain = abs(FTF);
phase = angle(FTF);
%%
figure('Color','w','Position',[100 100 1200 500]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% LEFT 
nexttile([2 1]);
stem(k*dt/tauR, h, 'k', 'filled', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('\tau / \tau_R');
ylabel('h');
title('Impulse response');
xlim([0 1.0]);
ylim([0 1.1]);
grid on;
box on;

%CENTER
nexttile([2 1]);
plot(real(FTF), imag(FTF), 'k', 'LineWidth', 1.5); 
hold on;
plot(real(FTF(540)), imag(FTF(540)), 'bo', 'MarkerFaceColor', 'k'); 
hold off;
xlabel('Re(F)');
ylabel('Im(F)');
title('Nyquist plot');
axis equal;
xlim([-2 2]);
ylim([-2 2]);
grid on;
box on;


%RIGHT
nexttile;
plot(omegaTauR, gain, 'k', 'LineWidth', 1.5);
hold on;
plot(omegaTauR(540), gain(540), 'bo', 'MarkerFaceColor', 'k'); 
hold off;
ylabel('|F|');
title('Bode plot');
xlim([0 25]);
ylim([0 2]);
grid on;
box on;

nexttile;
plot(omegaTauR, phase, 'k', 'LineWidth', 1.5);
hold on;
plot(omegaTauR(540), phase(540), 'bo', 'MarkerFaceColor', 'k'); 
hold off;
xlabel('\omega \tau_R');
ylabel('arg(F)');
xlim([0 25]);
ylim([-pi pi]);
yticks([-pi 0 pi]);
yticklabels({'-\pi','0','\pi'});
grid on;
box on;