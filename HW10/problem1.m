%% Alex Swenson ME573 HW10, Problem 1
function problem1()
clear; clc;

%% Problem Parameters

dx = 0.1; dy=dx;
% dt = 0.01;
x = -1:dx:1;
y = 0:dy:2;


% Numerical Parameters
kappa = 0.5;
nu=0.2;

% Meshgrid
[X, Y] = meshgrid(x,y);
X=X'; Y=Y';

% Exact Solution

[u_exact, v_exact] = get_exact(X,Y,kappa,nu);

subplot(2,1,1)
surf(X,Y,u_exact)
subplot(2,1,2)
surf(X,Y,v_exact)





end

function [u,v] = get_exact(X, Y, kappa, nu)

    x0 = 1;
    a0=0.001*kappa*exp((1+x0)*kappa); a1=a0; 

    u_num = -2*nu*(a1 + kappa*(exp(kappa*(X-x0))-exp(-kappa*(X-x0))*cos(kappa*Y)));
    v_num = 2*nu*(-kappa*(exp(-kappa*(X-x0)) + exp(kappa*(X-x0)))*cos(kappa*Y));
    
    
    denom = a0 + a1*X +(exp(-kappa*(X-x0)) + exp(kappa*(X-x0)))*cos(kappa*Y);
    
    u = u_num %/ denom;
    v = v_num %/ denom;
end