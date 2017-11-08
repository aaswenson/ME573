function problem1()
clear; clc; close all;
global dx nx ny f_exact f_guess f_source

% Set Domain and some important constants 
x_bound = 1; y_bound = 1;
dx = 0.05; dy = dx;
x = 0:dx:x_bound; y = 0:dy:y_bound;
nx = length(x); ny = length(y);
iteration_step = 20;
iteration_max = 800;

jac = zeros(1,(iteration_max / iteration_step));
gauss = zeros(1,(iteration_max / iteration_step));
N = 1:iteration_step:iteration_max;
% create mesh
[X, Y] = meshgrid(x,y);
X=X'; Y=Y';
% Get Analytic Solution
f_exact = X.*(1-X).*Y.*(1-Y);
% Set guess solution
f_guess = f_exact*0;
% get RHS vector
f_source = -2*X.*(1-X) - 2*Y.*(1-Y);

% call the comparisons
i = 1;
for n=1:iteration_step:iteration_max;
    [jac(i), gauss(i)] = run_comparison(n);
    i = i + 1;
end
% plot error results
figure
plot(N, jac,'g', N, gauss, 'r')
xlim([0,800])
ylim([0,0.07])
end
%% Run The Comparison Between Iteration Methods
function [jac, gauss] = run_comparison(N_iterations)
global dx nx ny f_exact f_guess f_source
f_save_jac = f_guess;
f_save_gauss = f_guess;
for it=1:N_iterations
    f_jac = jacobi_iteration(f_save_jac, f_source, nx, ny,dx);
    f_save_jac = f_jac;
    
    f_gauss = gauss_seidel(f_save_gauss, f_source, nx, ny,dx);
    f_save_gauss = f_gauss;
end
jac = max(max(abs(f_jac - f_exact)));
gauss = max(max(abs(f_gauss - f_exact)));
end
%% Jacobi Iteration Function
function f = jacobi_iteration(f, f_source, nx, ny,dx)


    f(2:nx-1,2:ny-1) = 0.25 * (f(1:nx-2,2:ny-1) + f(3:nx,2:ny-1)...
                             + f(2:nx-1,1:ny-2) + f(2:nx-1,3:ny))...
                             - 0.25 * (f_source(2:nx-1,2:ny-1)*dx^2);
end
%% Gauss-Seidel Iteration Function
function f = gauss_seidel(f, f_source, nx, ny, dx)
    for row=2:ny-1
       for col=2:nx-1
           f(row, col) = 0.25*(f(row, col+1) + f(row,col-1) + ...
               f(row+1,col) + f(row-1,col)) - 0.25*f_source(row,col)*dx^2;
       end
    end
end

