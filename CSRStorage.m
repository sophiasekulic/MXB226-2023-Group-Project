%% CSR Storage

[rb, c, v] = CSR(A);
%%
x0 = ones(33,1);
tol = 1e-10;
maxiters = 100000;

tic
[x_jacobi, converged_jacobi, k_jacobi, res_jacobi] = jacobitest(rb, v, c, B, x0, tol, maxiters)
toc
%%
k_jacobiplot = 0:k_jacobi;
figure
plot(k_jacobiplot,log(res_jacobi))
xlabel('Number of Iterations','Interpreter','latex')
ylabel('Residual','Interpreter','latex')
title(['Error of Jacobi Method at Each Iteration'],'Interpreter','latex')
%%
tic
[x_GS, converged_GS, k_GS, res_GS] = gaussseideltest(rb, v, c, B, x0, tol, maxiters)
toc
%%
k_GSplot = 0:k_GS;
figure
plot(k_GSplot,log(res_GS))
xlabel('Number of Iterations','Interpreter','latex')
ylabel('Residual','Interpreter','latex')
title(['Error of Gauss-Seidel Method at Each Iteration'],'Interpreter','latex')