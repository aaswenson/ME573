%% Alex Swenson ME573 HW10, Problem 2
function [x, f_ftcs, f_CN, f_analytic] = problem2()
clear; clc;

%% Problem Parameters
% physical constants
global kappa alpha

kappa = 0.02;

dx = 0.05;
dt = 0.05;
x = 0:dx:4;
L = 1.0;
N = length(x);
t_init = 0;
t_final = 2*(L*L/kappa);

alpha = kappa*dt / (dx*dx);

f_analytic = analytic(x,L,t_final);
f_ftcs = analytic(x,L,0);
f_CN = f_ftcs;

t = t_init;
while t < t_final;

    f_ftcs(2:N-1) = f_ftcs(2:N-1) - f_ftcs(2:N-1).*(dt/dx).*(f_ftcs(3:N) - f_ftcs(1:N-2))*0.5 + ...
               alpha*(f_ftcs(3:N) - 2*f_ftcs(2:N-1) + f_ftcs(1:N-2));
    f_ftcs(1) = analytic(x(1), L, t);
    f_ftcs(N) = analytic(x(N), L, t);
    t = t + dt;
end

Co = abs(f_analytic * (dt / dx));
Rec = Co / alpha;


end

function f = analytic(x,L,t)

    global kappa
    eterm = exp(-kappa*t/(L*L));
    f = -(2*kappa/L) * cosh(x/L) ./ (sinh(x/L) + eterm);
end

