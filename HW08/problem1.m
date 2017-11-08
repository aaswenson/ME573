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

i = 1;
for n=1:iteration_step:iteration_max;
    jac(i) = run_comparison(n);
    i = i + 1;
end

figure
plot(N, jac,'g')
xlim([0,800])
ylim([0,0.07])
end
function jac = run_comparison(N_iterations)
global dx nx ny f_exact f_guess f_source
f_save = f_guess;
for it=1:N_iterations
    f = jacobi_iteration(f_save, f_source, nx, ny,dx);
    f_save = f;
end
jac = max(max(abs(f - f_exact)));
end

function f = jacobi_iteration(f, f_source, nx, ny,dx)

    f(2:nx-1,2:ny-1) = 0.25 * (f(1:nx-2,2:ny-1) + f(3:nx,2:ny-1)...
                             + f(2:nx-1,1:ny-2) + f(2:nx-1,3:ny))...
                             - 0.25 * (f_source(2:nx-1,2:ny-1)*dx^2);
end

