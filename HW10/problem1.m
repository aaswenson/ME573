%% Alex Swenson ME573 HW10, Problem 1
function problem1()
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

f_init = analytic_solution(x,t_init,U,kappa,t_init);
f_exact = analytic_solution(x,t_final,U,kappa,t_init);

c = U*dt/dx;
alpha = kappa*dt/(dx*dx);

% Make the A matrix(2:N-1)
e = ones(N-1,1);
A = spdiags([(c+alpha)*e, (1-c-2*alpha)*e, alpha*e], -1:1, N-2, N-2);
f_save = f_init;
f_ftbs = f_init;

for i=1:dt:t_final
    for j=2:N-1
       f_ftbs(j) = (c+alpha)*f_save(j-1) + (1-c-alpha)*f_save(j) + alpha*f_save(j+1);
    end
    f_save = f_ftbs;
end

plot(x, f_ftbs, x, f_exact)

end

function A = eval_A(x, t, u, k)
     A = 0.5*erfc((x - u.*t) ./ sqrt(2*k*t)) + 0.5*exp(u*x/k).*erfc((x + u.*t)./sqrt(2*k*t));
end

function f = analytic_solution(x, t, u, k, t_0)

    C_i = pi;
    C_o = 2*pi;
    
    if t <= t_0
        f = C_i + (C_o - C_i)*eval_A(x,t,u,k); 
    else
        f = C_i + (C_o - C_i)*eval_A(x,t,u,k) - C_o*eval_A(x,t-t_0,u,k);
    end

end

