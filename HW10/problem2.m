%% Alex Swenson ME573 HW10, Problem 1
function problem2()
clear; clc;

%% Problem Parameters

dx = 0.05;
dt = 0.01;
x = 0:dx:4;
L = 1;
N = length(x);
t_init = 0.1;
t_final = 2.5;

% physical constants
U = 0.2;
kappa = 0.01;


f_exact = analytic_solution(x,t_final,U,kappa,t_init);


end

function f = analytic_solution(x, t, u, k, t_0)

    C_i = pi;
    C_o = 2*pi;
    
    if t > t_0
        f = C_i + (C_o - C_i)*eval_A(x,t,u,k) - C_o*eval_A(x,(t-t_0),u,k);
    else
        f = C_i + (C_o - C_i)*eval_A(x,t,u,k); 
    end

end

