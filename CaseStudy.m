%% Case Study
%% 
% 
%% Creating the Linear System

clear all
%%
%% GENERATING THE A MATRIX
a = (sqrt(2)+30)/15;
%% Imputting known values into a large matrix

value = -1 * [-2 1 1 ... %T1
    -3 1 1 1 ...    %T2 
    -3 1 1 1 ...    %T3
    -3 1 1 1 ...    %T4
    -3 1 1 1 ...    %T5
    -3 1 1 ...      %T6 
    -3 1 1 1 ...    %T7
    -4 1 1 1 1 ...  %T8
    -4 1 1 1 1 ...  %T9 
    -4 1 1 1 1 ...  %T10
    -4 1 1 1 1 ...  %T11
    -4 1 1 1 ...    %T12
    -3 1 1 1 ...    %T13
    -4 1 1 1 1 ...  %T14
    -4 1 1 1 1 ...  %T15
    -4 1 1 1 1 ...  %T16
    -4 1 1 1 1 ...  %T17
    -4 1 1 1 ...    %T18
    -3 1 1 1 ...    %T19
    -4 1 1 1 1 ...  %T20
    -4 1 1 1 1 ...  %T21
    -4 1 1 1 1 ...  %T22
    -4 1 1 1 1 ...  %T23
    -a 1 1 ...      %T24
    -3 1 1 1 ...    %T25
    -4 1 1 1 1 ...  %T26
    -4 1 1 1 1 ...  %T27
    -4 1 1 1 1 ...  %T28
    -a 1 1 ...      %T29
    -3 1 1 ...      %T30
    -4 1 1 1 ...    %T31
    -4 1 1 1 ...    %T32
    -a 1 1 ];       %T33
   
column_number = [1 2 7 ... %T1
    2 1 3 8 ...         %T2 
    3 2 4 9 ...         %T3
    4 3 5 10 ...        %T4
    5 4 6 11 ...        %T5
    6 5 12 ...          %T6 
    7 1 8 13 ...        %T7
    8 2 7 9 14 ...      %T8
    9 3 8 10 15 ...     %T9 
    10 4 9 11 16 ...    %T10
    11 5 10 12 17 ...   %T11
    12 6 11 18 ...      %T12
    13 7 14 19 ...      %T13
    14 8 13 15 20 ...   %T14
    15 9 14 16 21 ...   %T15
    16 10 15 17 22 ...  %T16
    17 11 16 18 23 ...  %T17
    18 12 17 24 ...     %T18
    19 13 20 25 ...     %T19
    20 14 19 21 26 ...  %T20
    21 15 20 22 27 ...  %T21
    22 16 21 23 28 ...  %T22
    23 17 22 24 29 ...  %T23
    24 18 23 ...        %T24
    25 19 26 30 ...     %T25
    26 20 25 27 31 ...  %T26
    27 21 26 28 32 ...  %T27
    28 22 27 29 33 ...  %T28
    29 23 28 ...        %T29
    30 25 31 ...        %T30
    31 26 30 32 ...     %T31
    32 27 31 33 ...     %T32
    33 28 32 ];         %T33

number_of_values_in_row = [3 4 4 4 4 3 4 5 5 5 5 4 4 5 5 5 5 4 4 5 5 5 5 3 4 5 5 5 3 3 4 4 3];

A = zeros(33,33);
tally = 1;
for i=1:33
    n = number_of_values_in_row(i);
    for j = tally:tally + n - 1
        c = column_number(j);
        A(i,c) = value(j);
    end
    tally = tally + n;
end

%% creating b matrix
b = zeros(33,1);
b([6 12 18]) = 40;
b([24 29 33]) = (4*sqrt(2))/3;
b(30:32) = 70;

%% solve for Temperature Values
Temperature = A\b
%% Storage Methods and Soutions
% Full Storage

%% Spy Plot of A
spy(A)
%% Factorising A
A_Full = cholesky_factorisation(A);

%% Spy Plot of Factorised A
spy(A_Full)
%% Solving for T
T_Full = General_Solve(A_Full, b)
%% Confirming Correct Solution
isequal(round(T_Full, 10), round(Temperature,10))
% Efficiency for Full Storage

%% Finding Bandwidths For A
kl = bandwidth(A_Full,'lower'); % Lower Bandwidth of A
ku = bandwidth(A_Full,'upper'); % Upper Bandwidth of A
kt = kl + ku + 1;         % Total Bandwidth
fprintf('Total Bandwidth of full storage: 13');
%% Calculating Fill in
A_Full_Fill = 250 - 141;
fprintf('Full Storage Fill in: 109');
%% Memory Usage
info_A_Full = whos('A_Full');
fprintf('Memory usage of Full Storage A: %.2f Bytes', info_A_Full.bytes);
%% Runtime for Solution
tic
General_Solve(A_Full, b);
toc
%% Number of Floating Point Operations


% Packed Storage

A_Packed = [];

n = size(A,1);
m = size(A,2);

for i = 1:n
    modified_A = A(i, 1:i);
    A_Packed = [A_Packed modified_A];
end

%% Apply Cholesky Factorisation for packed storage
A_Packed = Packed_cholesky_factorisation(A_Packed);

%% Solving for T
T_Packed = Packed_Solve(A_Packed, b)
%% Confirming Correct Solution
isequal(round(T_Packed, 10), round(Temperature,10))
% Efficiency for Packed Storage

%% Memory Usage
info_A_Packed = whos('A_Packed');
fprintf('Memory usage of Packed Storage A: %.2f Bytes', info_A_Packed.bytes);
%% Runtime
tic
Packed_Solve(A_Packed, b);
toc
%% Number of Floating Point Operations


% Band Storage

%% Performing Reverse Cuthill Mckee on matrix A 
band_perm = symrcm(A);

%% Reorder the matrix A and B using the permutation vector
Band_A_Perm_Full = A(band_perm, band_perm);
Band_B_Perm = b(band_perm);

%% Storing Only Upper Triangular Values
Band_A_Perm = triu(Band_A_Perm_Full);

%% Generating Spy Plots of Permutated A
spy(Band_A_Perm)
%% Factorising Permutated A
A_Band = Band_cholesky_factorisation(Band_A_Perm);

%% Generating Spy Plots of Permutated and Factorised A
spy(A_Band)
%% Calculating T with permutations using the A_Band
Band_T_Perm = General_Solve(A_Band', Band_B_Perm);

%% Reordering to obtain the original temperatures for each node
n = length(band_perm);
band_perm_inverse = zeros(1,n);

for i = 1:n
    band_perm_inverse(band_perm(i)) = i;
end

T_Band = Band_T_Perm(band_perm_inverse)
%% Confirming Correct Solution
isequal(round(T_Band, 10), round(Temperature,10))
% 
% Efficiency for Band Storage

%% Finding Bandwidths For Original Matrix A
ku_Band = bandwidth(A_Band,'upper'); % Lower Bandwidth of A
kt_Band = ku_Band + 1;         % Total Bandwidth
fprintf('Total Bandwidth of full storage: 7');

%% Calculating Fill in
A_Full_Fill = 168 - 87;
fprintf('Full Storage Fill in: 81');

%% Memory Usage
info_A_Band = whos('A_Band');
fprintf('Memory usage of Band Storage A: %.2f Bytes', info_A_Band.bytes);

%% Runtime for Solution
tic
General_Solve(A_Band, b);
toc

%% Number of Floating Point Operations
% Sparse Storage

%% Permutating Matrix A and vector B with AMD ordering
Sparse_Perm = symamd(A);
Sparse_A_Perm_Full = A(Sparse_Perm, Sparse_Perm);
Sparse_B_Perm = b(Sparse_Perm);

%% Converting to upper triangular matrix
Sparse_A_Perm = triu(Sparse_A_Perm_Full);

% Spy Plot of Permuted A
spy(Sparse_A_Perm)
%% Factorising Sparse version of A
Sparse_A = sparse_cholesky_factorisation(Sparse_A_Perm);

% Spy Plot of factorised and permutated A
spy(Sparse_A')
%% Calculating T
Sparse_T_Perm = General_Solve(Sparse_A, Sparse_B_Perm);

%% Reordering for original solution
zerop = 1:33;
invp(Sparse_Perm) = zerop;

T_Sparse = Sparse_T_Perm(invp)
%% Confirming Correct Solution
isequal(round(T_Sparse, 10), round(Temperature, 10))
% Efficiency for Sparse Storage

%% Finding Bandwidths For Original Matrix A
ku_Sparse = bandwidth(Sparse_A,'lower'); % Lower Bandwidth of A
kt_Sparse = ku_Sparse + 1;         % Total Bandwidth
fprintf('Total Bandwidth of Sparse storage: 30');
%% Calculating Fill in
A_Full_Fill = 140 - 87 %% Need to check second spy plot for fill in
fprintf('Sparse Storage Fill in: 53');
%% Memory Usage
info_A_Sparse = whos('Sparse_A');
fprintf('Memory usage of Sparse Storage A: %.2f Bytes', info_A_Sparse.bytes);
%% Runtime for Solution
tic
Sparse_A = sparse_cholesky_factorisation(Sparse_A_Perm);
General_Solve(Sparse_A, b);
toc
%% Number of Floating Point Operations

% CSR Storage

%% Converting A into Compressed Sparse Row Format
[rb, c, v] = CSR(A)
%% Iterative Methods
% Jacobi

x0 = ones(33,1); % Initial guess
tol = 1e-14; % Tolerance
maxiters = 10000; % Number of maximum iterations

% Jacobi Method for CSR Storage
[T_Jacobi, Converged_Jacobi, k_Jacobi, Res_Jacobi] = Jacobi_CSR(rb, v, c, b, x0, tol, maxiters)
%% Confirming Correct Solution
isequal(round(T_Jacobi, 10), round(Temperature, 10))
%% 
% *Efficiency for Jacobi Method using CSR Storage*

%% Runtime for Solution
tic
[T_Jacobi, Converged_Jacobi, k_Jacobi, Res_Jacobi] = Jacobi_CSR(rb, v, c, b, x0, tol, maxiters);
toc
%% Number of Iterations
k_Jacobi
%% Number of Floating Point Operations

%% 
% *Tolerance for Jacobi*

% Finding the tolerance required to produce a solution visually indistiguishable from that 
% produced by the direct methods. Assume 'visually indistiguishable' means to 4 decimal places. 

x0 = ones(33,1); % Initial guess
maxiters = 10000; % Number of maximum iterations

tol1 = 1e-7;
[T_JTol1, ~, ~, ~] = Jacobi_CSR(rb, v, c, b, x0, tol1, maxiters);
isequal(round(T_JTol1, 4), round(Temperature, 4))
tol2 = 1e-8;
[T_JTol2, ~, ~, ~] = Jacobi_CSR(rb, v, c, b, x0, tol2, maxiters);
isequal(round(T_JTol2, 4), round(Temperature, 4))
tol3 = 1e-9;
[T_JTol3, ~, ~, ~] = Jacobi_CSR(rb, v, c, b, x0, tol3, maxiters);
isequal(round(T_JTol3, 4), round(Temperature, 4))
% Gauss-Seidel

x0 = ones(33,1); % Initial guess
tol = 1e-14; % Tolerance
maxiters = 10000; % Number of maximum iterations

% Gauss-Seidel Method for CSR Storage
[T_GS, Converged_GS, k_GS, Res_GS] = GaussSeidel_CSR(rb, v, c, b, x0, tol, maxiters)
%% Confirming Correct Solution
isequal(round(T_GS, 10), round(Temperature, 10))
%% 
% *Efficiency for Gauss-Seidel Method using CSR Storage*

%% Runtime for Solution
tic
[T_GS, Converged_GS, k_GS, Res_GS] = GaussSeidel_CSR(rb, v, c, b, x0, tol, maxiters);
toc
%% Number of Iterations
k_GS
%% Number of Floating Point Operations

%% 
% *Tolerance for Gauss-Seidel*

% Finding the tolerance required to produce a solution visually indistiguishable from that 
% produced by the direct methods. Assume 'visually indistiguishable' means to 4 decimal places.

x0 = ones(33,1); % Initial guess
maxiters = 10000; % Number of maximum iterations

tol1 = 1e-7;
[T_GSTol1, ~, ~, ~] = GaussSeidel_CSR(rb, v, c, b, x0, tol1, maxiters);
isequal(round(T_GSTol1, 4), round(Temperature, 4))
tol2 = 1e-8;
[T_GSTol2, ~, ~, ~] = GaussSeidel_CSR(rb, v, c, b, x0, tol2, maxiters);
isequal(round(T_GSTol2, 4), round(Temperature, 4))
tol3 = 1e-9;
[T_GSTol3, ~, ~, ~] = GaussSeidel_CSR(rb, v, c, b, x0, tol3, maxiters);
isequal(round(T_GSTol3, 4), round(Temperature, 4))
% SOR

x0 = ones(33,1); % Initial guess
tol = 1e-14; % Tolerance
maxiters = 10000; % Number of maximum iterations

% Calculating omega

D = diag(diag(A)); L = tril(A,-1); U = triu(A,1);
T_J = -inv(D)*(L+U); % Jacobi iteration matrix
rhoT_J = max(abs(eig(T_J))); % spectral radius of T_J
omega = 2 / (1+sqrt(1-(rhoT_J)^2)) % calculating optimal omega using Laplace's equation
% SOR Method for CSR Storage
[T_SOR, Converged_SOR, k_SOR, Res_SOR] = SOR_CSR(rb, v, c, b, x0, tol, maxiters, omega)
%% Confirming Correct Solution
isequal(round(T_SOR, 10), round(Temperature, 10))
%% 
% *Efficiency for SOR Method using CSR Storage*

%% Runtime for Solution
tic
[T_SOR, Converged_SOR, k_SOR, Res_SOR] = SOR_CSR(rb, v, c, b, x0, tol, maxiters, omega);
toc
%% Number of Iterations
k_SOR
%% Number of Floating Point Operations

%% 
% *Tolerance for SOR*

% Finding the tolerance required to produce a solution visually indistiguishable from that 
% produced by the direct methods. Assume 'visually indistiguishable' means to 4 decimal places.

x0 = ones(33,1); % Initial guess
maxiters = 10000; % Number of maximum iterations

tol1 = 1e-7;
[T_SORTol1, ~, ~, ~] = SOR_CSR(rb, v, c, b, x0, tol1, maxiters, omega);
isequal(round(T_SORTol1, 4), round(Temperature, 4))
tol2 = 1e-8;
[T_SORTol2, ~, ~, ~] = SOR_CSR(rb, v, c, b, x0, tol2, maxiters, omega);
isequal(round(T_SORTol2, 4), round(Temperature, 4))
%% 
% *Effect of omega on rate of convergence*

% Calculating the spectral radius with varying omega

D = diag(diag(A)); L = tril(A,-1); U = triu(A,1);
N = 1000;
omega = linspace(0.01,2,N);
rho_TSOR = zeros(N,1);
for i = 1:N
    TSOR = inv(D/omega(i) + L)*(((1-omega(i))/omega(i))* D - U);
    rho_TSOR(i) = max(abs(eig(TSOR)));
end

T_J = -inv(D)*(L+U); % Jacobi iteration matrix
rhoT_J = max(abs(eig(T_J))); % spectral radius of T_J
omega_opt = 2 / (1+sqrt(1-(rhoT_J)^2)); % calculating optimal omega using Laplace's equation
TSOR_opt = inv(D/omega_opt + L)*(((1-omega_opt)/omega_opt)* D - U);
rho_opt = max(abs(eig(TSOR_opt)));

figure
plot(omega, rho_TSOR, 'DisplayName', '')
hold on
plot(omega_opt, rho_opt, 'r.', 'MarkerSize', 20, 'DisplayName','Optimal Omega')
xlabel('Omega','Interpreter','latex')
ylabel('rho(TSOR)','Interpreter','latex')
title(['Spectral Radius of SOR Method with Varying Omega'],'Interpreter','latex')
legend
hold off
% Calculating the number of iterations required with varying omega

tol = 1e-14;
maxiters = 10000;
n = 1500;
omega2 = linspace(0.45,1.95,n);
k_omega = zeros(size(omega2));

for i = 1:n
    [~, ~, k_omega(i), ~] = SOR_CSR(rb, v, c, b, x0, tol, maxiters, omega2(i));
end

figure
plot(omega2, k_omega, 'DisplayName','')
hold on
plot(omega_opt, k_SOR, 'r.', 'MarkerSize', 20, 'DisplayName','Optimal Omega')
xlabel('Omega','Interpreter','latex')
ylabel('Number of Iterations','Interpreter','latex')
title(['Number of Iterations Required to Converge with Varying Omega'],'Interpreter','latex')
legend
hold off
% Conjugate Gradient

x0 = ones(33,1); % Initial guess
tol = 1e-14; % Tolerance
maxiters = 10000; % Number of maximum iterations

% Conjugate Gradient Method for CSR Storage
[T_CG, Converged_CG, k_CG, Res_CG] = CG_CSR(rb, v, c, b, x0, tol, maxiters)
%% Confirming Correct Solution
isequal(round(T_CG, 10), round(Temperature, 10))
%% 
% *Efficiency for Conjugate Gradient Method using CSR Storage*

%% Runtime for Solution
tic
[T_CG, Converged_CG, k_CG, Res_CG] = CG_CSR(rb, v, c, b, x0, tol, maxiters);
toc
%% Number of Iterations
k_CG
%% Number of Floating Point Operations

%% 
% *Tolerance for Conjugate Gradient*

% Finding the tolerance required to produce a solution visually indistiguishable from that 
% produced by the direct methods. Assume 'visually indistiguishable' means to 4 decimal places.

x0 = ones(33,1); % Initial guess
maxiters = 10000; % Number of maximum iterations

tol1 = 1e-6;
[T_CGTol1, ~, ~, ~] = CG_CSR(rb, v, c, b, x0, tol1, maxiters);
isequal(round(T_CGTol1, 4), round(Temperature, 4))
tol2 = 1e-7;
[T_CGTol2, ~, ~, ~] = CG_CSR(rb, v, c, b, x0, tol2, maxiters);
isequal(round(T_CGTol2, 4), round(Temperature, 4))
%% 
% *Comparison of the Iterative Methods*

k_JacobiPlot = 0:k_Jacobi;
k_GSPlot = 0:k_GS;
k_SORPlot = 0:k_SOR;
k_CGPlot = 0:k_CG;
figure
plot(k_JacobiPlot,log(Res_Jacobi),'r-','DisplayName','Jacobi')
hold on
plot(k_GSPlot,log(Res_GS),'b-','DisplayName','Gauss-Seidel')
plot(k_SORPlot,log(Res_SOR),'m-','DisplayName','SOR')
plot(k_CGPlot,log(Res_CG),'g-','DisplayName','Conjugate Gradient')
xlabel('Number of Iterations','Interpreter','latex')
ylabel('Log of the Residual Norm','Interpreter','latex')
title(['Residual Norm of the Iterative Methods'],'Interpreter','latex')
legend
hold off
%% Visualisation

%Intialising electronic component matrix
T = zeros(7,7);

%Assigning temperatures to specifc nodes
T(1, 1:4) = 70;
T(4:7, 7) = 40;
T(7, 1:6) = Temperature(1:6);
T(6, 1:6) = Temperature(7:12);
T(5, 1:6) = Temperature(13:18);
T(4, 1:6) = Temperature(19:24);
T(3, 1:5) = Temperature(25:29);
T(2, 1:4) = Temperature(30:33);
T(T == 0) = NaN;

%Generating the Mesh Grid for component
[X, Y] = meshgrid(0:0.01:0.06, 0.06:-0.01:0);

Temperature = T(~isnan(T));
x = X(~isnan(T));
y = Y(~isnan(T));

%Triangulating to connect vertices
tri = delaunay(x, y);

%Visualising the mesh with numbered triangles
figure
trimesh(tri, x, y)
axis image

%Removing unused triangles
tri([29, 44, 36],:) = [];

%Numbering the triagnles
for i = 1:size(tri,1);
    indices = tri(i,:);
    vertices = [x(indices) y(indices)];
    centroid = mean(vertices);
    text(centroid(1), centroid(2), num2str(i));
end

%Creating the final visualisation
figure
trisurf(tri, x, y, Temperature)
view(2)
shading interp
colormap cool
colorbar
axis image
grid off
%% Effect of Ambient Temperature

%% Point (0.03, 0.03) corresponds to node 22

%% Constructing new B vector where the ambient temperature is kept as an unknown variable
%% 'a' corresponds to the ambient temperature
syms a
BAmbient = zeros(33,1);
BAmbient = symmatrix(BAmbient);
BAmbient(6) = 40;
BAmbient(12) = 40;
BAmbient(18) = 40;
BAmbient(24) = a.*sqrt(2)/15;
BAmbient(29) = a.*sqrt(2)/15;
BAmbient(30) = 70;
BAmbient(31) = 70;
BAmbient(32) = 70;
BAmbient(33) = a.*sqrt(2)/15;
BAmbient = symmatrix2sym(BAmbient);
sympref('FloatingPointOutput',true);
Tsym = A \ BAmbient
Node22 = Tsym(22) % Node corresponding to T(0.03, 0.03)
%% Calculating the ambient temperature which makes T(0.03, 0.03) 50 degrees
LowTemp = Node22 == 50
LowestTemp = solve(LowTemp, a)
%% Calculating the ambient temperature which makes T(0.03, 0.03) 55 degrees
HighTemp = Node22 == 55
HighestTemp = solve(HighTemp, a)
% Check by calculating the temperature at node 22 with ambient temps ranging from
% 0 to 60 degrees

AmbientTemps = 0:1:60;
Node22Temps = zeros(1,length(AmbientTemps));

for i = 0:length(AmbientTemps)-1
    B1 = BAmbient;
    B1 = subs(B1, a, i);
    CheckNode22 = A \ B1;
    Node22Temps(i+1) = CheckNode22(22);
end

AmbientTemps, Node22Temps
index = find(Node22Temps < 50 | Node22Temps > 55)
OutsideRangeTemps = AmbientTemps(index)
%% 
%