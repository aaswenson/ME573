%% Alex Swenson ME573 HW10, Problem 2
function problem2()
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

%% Build CN matrices and run solvers


% Make the A matrix
A = update_A(f_CN, N, dx, dt);
% Make the B matrix
b = update_b(f_CN, N, dt, dx);

for i=t_init:dt:t_final
    f_ftcs(2:N-1) = f_ftcs(2:N-1) - f_ftcs(2:N-1).*(dt/dx).*(f_ftcs(3:N) - f_ftcs(1:N-2))*0.5 + ...
               alpha*(f_ftcs(3:N) - 2*f_ftcs(2:N-1) + f_ftcs(1:N-2));
    f_CN(2:N-1) = thomas(A,b);
    A = update_A(f_CN, N, dx, dt);
    b = update_b(f_CN, N, dt, dx);
end

Co = abs(f_analytic * (dt / dx));
Rec = Co / alpha;
% figure(1)
% semilogy(x, Rec, x, Co)
% figure(2)
figure(1)
plot(x,f_ftcs, x, f_analytic)
figure(2)
plot(x,f_CN)


end

function f = analytic(x,L,t)

    global kappa
    eterm = exp(-kappa*t/(L*L));
    f = -(2*kappa/L) * cosh(x/L) ./ (sinh(x/L) + eterm);
end

function b = update_b(f, N, dt, dx)
    global alpha
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    b = (f(2:N-1)*(dt/dx) + 2*alpha).*f(1:N-2) + 4*(1-alpha)*f(2:N-1) + (2*alpha - f(2:N-1)*(dt/dx)).*f(3:N);

end

function A = update_A(f, N, dx, dt)

global alpha 
% Make the A matrix(2:N-1)
e = ones(N,1);
A = spdiags([-(f'*(dt/dx) + 2*alpha).*e, 4*(1+alpha)*e, (2*alpha - f'*(dt/dx) ).*e], -1:1, N-2, N-2);

end
