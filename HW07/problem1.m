%% Alex Swenson HW7 (Computer Project 5c)
function problem1()
clear; clc; close all;
% Set Domain and some important constants 
x_bound = 1; y_bound = 1;
dx = 0.05; dy = dx;
x = 0:dx:x_bound; y = 0:dy:y_bound;
global k 
k = 0.1;
nx = length(x); ny = length(y);
dx2 = dx*dx; dy2 = dy*dy;
dt = (dx2 / (2*k));
t_end = 40 * dt;
nt = t_end / dt;
alphax = k*dt / dx2; alphay = k*dt / dy2;

% create mesh
[X, Y] = meshgrid(x,y);
X=X'; Y=Y';
% Set Initial Conditions
f = X .* (1-X.^5) .* Y .* (1-Y);
f(1,:) = 0;
f(:,1) = 0;

%% Get Exact Solution
Ntrunc = 50;
f_exact = get_exact(X, Y, Ntrunc, t_end, nx, ny);
f_exact(1,:) = 0;
f_exact(:,1) = 0;

% Part A) L_inf Error Between f_init, f_exact @ 0
A_exact = get_exact(X, Y, Ntrunc, 0, nx, ny);
L_inf = max(max(abs(f - A_exact)));
out = ['The L_inf error between fexact(t=0) and f_initial is: ', num2str(L_inf)];
disp(out);

%% Solve the Problem Using ADI

% Make The A - matrices
ex = ones(nx-2, 1);
ey = ones(ny-2,1);
A_x = spdiags([-alphax*ex 2*(1+alphax)*ex -alphax*ex], -1:1, nx-2, nx-2);
A_y = spdiags([-alphay*ey 2*(1+alphay)*ey -alphay*ey], -1:1, ny-2, ny-2);


% Allocate the f, b matrices
f_inter = f;
f_save = f;

for t=1:2
    step = mod(t, 2);

    if step == 1
        for j=2:nx-1
            b = (2 + alphay*dy2)*f(:,j);
            f_inter(2:ny-1,j) = A_x \ b(2:ny-1);
        end
        for i=2:ny-1
            b = (2 + alphax*dx2)*f_inter(i,:);
            f_save(i,2:nx-1) = A_y \ b(2:nx-1)';
        end
    elseif step == 0
        for i=2:ny-1
            b = (2 + alphax*dx2)*f(i,:);
            f_inter(i,2:nx-1) = A_y \ b(2:nx-1)';
        end
        for j=2:nx-1
            b = (2 + alphay*dy2)*f_inter(:,j);
            f_save(2:ny-1,j) = A_x \ b(2:ny-1);
        end
    end
    f = f_save;
end
%% Plot Results

tend = dt * nt;

%% Plotting Results
figure(1)
surf(X,Y,f_exact)
xlabel('X [-]')
ylabel('Y [-]')
title(['Exact Solution at T= ', num2str(tend)])
saveas(gcf, './writeup/exact_solution.png')

figure(2)
surf(X,Y,f);
xlabel('X [-]')
ylabel('Y [-]')
title(['ADI Solution at T= ', num2str(tend)])
saveas(gcf, './writeup/adi_solution.png')

figure(3)
surf(X,Y,abs(f - f_exact));
xlabel('X [-]')
ylabel('Y [-]')
title(['Solution Error (ADI vs. Exact) at T= ', num2str(tend)])
saveas(gcf, './writeup/solution_error.png')

end % End of Main Program


% This function plots the exact solution to the 2D diffusion PDE
function exactf = get_exact(X, Y, trunc, t_end, nx, ny)
global k
exactf = zeros(length(nx), length(ny));
for n=1:trunc
    for m=1:trunc
    exactf=exactf-120*(-n^4*pi^4*(-1)^n + 12*n^2*pi^2*(-1)^n + 24 + 24*(-1)^(1+n))...
        * (-2 + 2*(-1)^m) / (n^7*pi^10*m^3) ...
        * sin(n*pi*X) .*sin(m*pi*Y)*exp(-(n^2+m^2)*pi^2*k*t_end);
    end
end
end
