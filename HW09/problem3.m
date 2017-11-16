%% Numerical Solutions to the Convection equation
%
%       df/dt + U*(df/dx) = 0
%
function problem3()
clear; clc;
%% Set constants, problem space

U = pi;
dt = 0.01; tf = 0.5;
nt = (tf / dt);
dx = 0.5;
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


%% Solve using CN Method
f = f_init;

e = ones(N-1,1);
A = spdiags([-C_0*e, 4*e, C_0*e], -1:1, N-2, N-2);
b = C_0*f(1:N-2) + 4*f(2:N-1) - C_0*f(3:N);
b(1) = b(1) + 2*C_0;
b(N-2) = b(N-2) + 2*C_0;

for i=1:nt
   f(2:N-1) = TDMAsolver(A,b);
   b = C_0*f(1:N-2) + 4*f(2:N-1) - C_0*f(3:N);
   b(1) = b(1) + 2*C_0;
   b(N-2) = b(N-2) + 2*C_0;
end

plot(x,f, x, f_analytic)
end

function f= exact(x,t,U)
    %% get exact solution f(x,t)
    f = (erf((1-(x-U*t))/0.25) - erf((1+(x-U*t))/0.25));
end