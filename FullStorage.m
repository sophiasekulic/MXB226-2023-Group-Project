clear all
%%
% Full Storage
tic
G_full = cholesky_factorisation(A)
isequal(round(A, 12), round(G_full*G_full', 12))
z_full = forward_substitution(G_full, B);
x_full = backward_substitution(G_full', z_full)
toc
isequal(round(A*x_full, 10), round(B,10))
%whos