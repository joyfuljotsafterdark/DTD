clear all;
close all;
clc;

%% ========== load data ============ %%

% raw_data = readmatrix("data\test3_FTF_output.csv");
% 
% freq = raw_data(3:end,2);

%% ========== load data ============ %%
N = 5000;
% omega = 0:N-1;
freq = 0:N-1;
omega = 2*pi*freq;

k = 0:N-1;

dt = 1;

A = zeros(N,N);
%%
for m = 1:N
    for n = 1:N
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


% offdiag_col = Gcol - diag(diag(Gcol));
% offdiag_row = Grow - diag(diag(Grow));
% 
% err_col = norm(offdiag_col, 'fro');
% err_row = norm(offdiag_row, 'fro');
% 
% disp(offdiag_col)
% disp(offdiag_row)
% disp(err_col)
% disp(err_row)
% disp(is_col_orthogonal)
% disp(is_row_orthogonal)