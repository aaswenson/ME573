%% Numerical Solutions to the Convection equation
%
%       df/dt + U*(df/dx) = 0
function problem2()
clear; clc;
%% Set constants, problem space

U = pi;
dt = 0.01; tf = 0.5;
nt = (tf / dt);
dx = 0.05;
x = -5:dx:5;
N = length(x) - 1;

% Initial condition
f_init = exact(x,0,U);
% Final analytic solution
f_analytic = exact(x,tf,U);
% Check for numerical stability
C_0 = U * dt / dx;
if C_0 > 1
    disp(['Warning, C_0 = ', num2str(C_0, 3), ' your solution is numerically unstable!']);
else 
    disp(['C_0 = ', num2str(C_0,3)]);
end

%% Make the Coefficient Matrices

e = ones(N-2,1);
% A_ftcs = spdiags([0.5*C_0*e, e, -0.5*C_0*e], -1:1,N-2,N-2);
A_ftbs = spdiags([e*C_0, e*(1-C_0), e*0], -1:1,N-2,N-2);

full(A_ftbs)

%% Solve using FTBS method
f_ftbs = f_init;
for i=1:nt
    f_save(2:N-1) = f_ftbs(2:N-1) * A_ftbs;
    f_ftbs(2:N-1) = f_save(2:N-1);
end

plot(x,f_ftbs, x, f_analytic)

end

function f= exact(x,t,U)
    %% get exact solution f(x,t)
    f = (erf((1-(x-U*t))/0.25) - erf((1+(x-U*t))/0.25));
end