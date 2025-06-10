%% Sparse Storage

tic
U = A;
%%
p = symamd(U)
%%
U = U(p,p); %permutation
%%
isequal(U, U')
%%
B_sparse = B(p) % permuting B vector
%%
% Only upper portion needs to be stored ?
U2 = triu(U)
%%
U_sparse = sparse(U2)
%%
check_sparse = chol(U_sparse); % MATLAB built in function to check
%%
isequal(round(check_sparse'*check_sparse,10), round(sparse(U),10))
check_sparse2 = check_sparse';
%%
%G2 = cholesky_factorisation_sparse(U2)
G_sparse = cholesky_factorisation_sparse2(U_sparse);
%%
% Checking it is the same as built in function
isequal(round(G_sparse, 10), round(check_sparse2,10))
%%
% Checking it is the correct decomposition GG'=A
isequal(round(G_sparse*G_sparse',10), round(sparse(U),10))
%%
z_sparse = forward_substitution(G_sparse, B_sparse);
x_sparse = backward_substitution(G_sparse', z_sparse)
toc
%%
isequal(round(U*x_sparse, 10), round(B_sparse, 10))
%%
zerop = 1:33;
%p2 = linspace[1,33,33];
%invp(p)=zerop(1:33)
invp(p) = zerop % inverse of permutation vector to restore order of temps by node
% reordering B vector to original ordering
x_sparse2 = x_sparse(invp)
isequal(round(x_sparse2, 10), round(x_full, 10))
%%
figure, spy(U2')
figure, spy(G_sparse)
figure, spy(tril(A))
figure, spy(G_full)