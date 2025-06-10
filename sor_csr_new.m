function [x, converged, k] = sor_csr_new(rb, c, v, b, x0, tol, maxiters, omega)

% Initialize
n = length(rb) - 1; % Number of rows/columns
k = 0;
x = x0;
normb = norm(b);

% Initialize Ax using CSR format
Ax = zeros(n, 1);
for i = 1:n
    Ax(i) = v(rb(i):rb(i+1)-1) * x(c(rb(i):rb(i+1)-1));
end

% Compute the initial residual
res = norm(b - Ax) / normb;

%% Perform the iteration until res <= tol or maximum iterations taken

while k < maxiters && res > tol
    for i = 1:n
        x(i) = b(i) + (1 - omega) * x(i) * v(rb(i):rb(i+1)-1) * x(c(rb(i):rb(i+1)-1)) / omega;
        for j = rb(i):(rb(i+1)-1)
            if c(j) ~= i
                x(i) = x(i) - v(j) * x(c(j));
            end
        end
        x(i) = omega * x(i) / v(rb(i):rb(i+1)-1) * x(c(rb(i):rb(i+1)-1));
    end

    % Update Ax using CSR format
    for i = 1:n
        Ax(i) = v(rb(i):rb(i+1)-1)' * x(c(rb(i):rb(i+1)-1));
    end

    % Compute the residual
    res = norm(b - Ax) / normb;
    
    k = k + 1;
end

%% Check convergence

if res <= tol
    converged = true;
else
    converged = false;
end
end