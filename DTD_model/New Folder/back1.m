clear all;
close all;
clc;

%% ========== load data ============ %%

raw_data = readmatrix("data\test3_FTF_output.csv");

freq = raw_data(3:end,2);

%% ========== load data ============ %%
omega = 2*pi*freq;

k = 0:14;

dt = 5e-5;

A = zeros(15,15);
%%
for m = 1:15
    for n = 1:15
        A(m,n) = exp(-1*i * omega(m)*dt*k(n));
    end
end
%% check orthnormal%%

Gcol = A' * A; 
Grow = A * A';

tol = 1e-10;

is_col_orthogonal = norm(Gcol - diag(diag(Gcol)), 'fro') < tol;
is_row_orthogonal = norm(Grow - diag(diag(Grow)), 'fro') < tol;

disp(is_col_orthogonal)
disp(is_row_orthogonal)

%%
kappa = cond(A,2);

%%