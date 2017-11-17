%% Numerical Solutions to the Convection equation
%
%       df/dt + U*(df/dx) = 0
%
function [x, f_analytic, f] = problem3(dt, dx, tf)
%% Set constants, problem space

U = pi;
nt = (tf / dt);
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

% save initial boundary conditions
bc1 = f_init(1);
bc2 = f_init(N-1);


%% Solve using CN Method

% Make the A matrix(2:N-1)
e = ones(N-1,1);
A = spdiags([-C_0*e, 4*e, C_0*e], -1:1, N-2, N-2);

% make the first RHS matrix
b = update_b(f_init, C_0, N, bc1, bc2);
f = f_init;
for i=1:nt
   f(2:N-1) = A\b';
   b = update_b(f, C_0, N, bc1, bc2);
end

end

function b = update_b(f, C_0, N, bc1, bc2)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    b = C_0*f(1:N-2) + 4*f(2:N-1) - C_0*f(3:N);
    b(1) = -C_0*f(3) + 4*f(2) + C_0*(bc1 + bc1);
    b(N-2) = -C_0*(bc2 + bc2) + 4*f(N-2) + C_0*f(N-3);
end

function f= exact(x,t,U)
    %% get exact solution f(x,t)
    f = (erf((1-(x-U*t))/0.25) - erf((1+(x-U*t))/0.25));
end