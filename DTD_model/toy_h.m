clear all;
close all;
clc;

%%

L = 5;


v1 = (1:L)';        
v2 = 1:L;          


M = v1 * v2;
disp(M);

%%
new_matrix = zeros(L,L);

for i = 1:L
        new_matrix(i,i:L) = M(i,1:(L-i+1));
end

disp(new_matrix);

%%
r = sum(new_matrix,1);
disp(r)

%%
v1_restored = zeros(L, 1);

for j = 1:L
    
    sum_known = 0;
    for i = 1:(j-1)
        sum_known = sum_known + v1_restored(i) * v2(j - i + 1);
    end
    
    v1_restored(j) = (r(j) - sum_known) / v2(1);
end

disp(v1_restored');








%%