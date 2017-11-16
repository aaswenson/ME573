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


%% Solve using FTBS and FTCS method
f_ftbs = f_init;
f_ftcs = f_init;
f_save_forward = f_ftbs;
f_save_central = f_ftcs;
for i=1:nt
    for j=2:N-1
       f_ftbs(j) = f_save_forward(j-1)*C_0 + (1-C_0)*f_save_forward(j);
       f_ftcs(j) = (C_0*0.5)*f_save_central(j-1) + f_save_central(j) - (C_0*0.5)*f_save_central(j+1);
    end
    f_save_forward = f_ftbs;
    f_save_central = f_ftcs;
end
figure(1)
plot(x, f_analytic, 'r', x, f_ftcs, 'g', x, f_ftbs, 'b')
title(['Comparison of Linear Advection Schemes'])
xlabel(['X [-]'])
ylabel(['f(X) [-]'])
legend('Analytic', 'FTCS', 'FTBS')
saveas(gcf,['./writeup/p2_compare_schemes.png'])

end

function f= exact(x,t,U)
    %% get exact solution f(x,t)
    f = (erf((1-(x-U*t))/0.25) - erf((1+(x-U*t))/0.25));
end