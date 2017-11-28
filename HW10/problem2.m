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
    f_ftcs = ftcs(x, t, f_ftcs, N, dt, dx, L);
    f_CN = CN(x, t, f_CN, N, dt, dx, L);
    t = t + dt;
end

end

function f = analytic(x,L,t)

    global kappa
    eterm = exp(-kappa*t/(L*L));
    f = -(2*kappa/L) * cosh(x/L) ./ (sinh(x/L) + eterm);
end

function f_ftcs = ftcs(x, t, f_ftcs, N, dt, dx, L)
    % execute one time step with FTCS scheme
    global alpha
    f_ftcs(2:N-1) = f_ftcs(2:N-1) - f_ftcs(2:N-1).*(dt/dx).*(f_ftcs(3:N) - f_ftcs(1:N-2))*0.5 + ...
               alpha*(f_ftcs(3:N) - 2*f_ftcs(2:N-1) + f_ftcs(1:N-2));
    f_ftcs(1) = analytic(x(1), L, t);
    f_ftcs(N) = analytic(x(N), L, t);

end

function f_CN = CN(x, t, f_CN, N, dt, dx, L)
    % Execute one time step with CN scheme
    A = update_A(f_CN, N, dt, dx);
    b = update_b(f_CN, N, dt, dx);
    f_CN(2:N-1) = thomas(A,b(2:N-1));
    f_CN(1) = analytic(x(1), L, t);
    f_CN(N) = analytic(x(N), L, t);
    
end


function A = update_A(f, N, dt, dx)
    global alpha
    e = ones(N-1, 1);
    % set ghost point to make matrix
    f(N+1) = 1;
    A = spdiags([-((f(3:N+1)'*dt/dx) +2*alpha).*e, 4*(1+alpha)*e, (f(1:N-1)'*(dt/dx) - 2*alpha).*e], -1:1, N-2, N-2);
end

function b = update_b(f, N, dt, dx)
    global alpha
    b = ones(1,N);
    b(2:N-1) = (f(2:N-1) * (dt/dx) + 2*alpha).*f(1:N-2) + ...
                4*(1-alpha)*f(2:N-1) + ...
               (2*alpha - f(2:N-1)*(dt/dx)).*f(3:N);
    b(2) = b(2) + (f(2)*dt / dx + 2*alpha)*f(1);
    b(N-1) = b(N-1) + (2*alpha - f(N-1)*dt / dx)*f(N);
    
end


