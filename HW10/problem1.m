%% Alex Swenson ME573 HW10, Problem 1
function [x, f_ftbs, f_ftcs, f_CN, f_exact] = problem1()
clear; clc;

%% Problem Parameters

dx = 0.05;
dt = 0.01;
x = -5:dx:5;
N = length(x);
t_init = 0.1;
t_final = 2.5;

% physical constants
U = 0.2;
kappa = 0.01;

f_init = analytic_solution(x,t_init,U,kappa);
f_exact = analytic_solution(x,t_final,U,kappa);

c = U*dt/dx;
alpha = kappa*dt/(dx*dx);

f_ftbs = f_init;
f_ftcs = f_init;
f_CN = f_init;

bc1 = f_init(1);
bc2 = f_init(N-1);
%% Set Up CN matrices

% Make the A matrix(2:N-1)
e = ones(N-1,1);
A = spdiags([-(c+2*alpha)*e, 4*(1+alpha)*e, (c-2*alpha)*e], -1:1, N-2, N-2);

% Make the B matrix
b = update_b(f_CN,c,N,alpha, bc1, bc2);


%% March Forward in Time
for i=1:dt:t_final
    f_ftbs(2:N-1) = (c+alpha)*f_ftbs(1:N-2) + (1-c-2*alpha)*f_ftbs(2:N-1) + alpha*f_ftbs(3:N);
    f_ftcs(2:N-1) = (0.5*c + alpha)*f_ftcs(1:N-2) + (1 - 2*alpha)*f_ftcs(2:N-1) - (0.5*c - alpha)*f_ftcs(3:N);
    f_CN(2:N-1) = thomas(A,b);
    b = update_b(f_CN, c, N, alpha, bc1, bc2);
    
end

end

function A = eval_A(x, t, u, k)
     A = 0.5*erfc((x - u*t) ./ sqrt(2*k*t)) + 0.5*exp(u*x/k).*erfc((x + u*t)./sqrt(2*k*t));
end

function f = analytic_solution(x, t, u, k)
    C_i = pi;
    C_o = 2*pi;
    f = C_i + (C_o - C_i)*eval_A(x,t,u,k);
end

function b = update_b(f, c, N, alpha, bc1, bc2)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    b = (c+2*alpha)*f(1:N-2) + 4*(1-alpha)*f(2:N-1) - (c-2*alpha)*f(3:N);
    b(1) = b(1) + (c+2*alpha)*bc1;
    b(N-2) = b(N-2) - (c-2*alpha)*bc2;

end
