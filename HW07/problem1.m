%% Alex Swenson HW7 (Computer Project 5c)
function problem1()
% Set Domain and some important constants 
x_bound = 1; y_bound = 1;
dx = 0.05; dy = dx;
global k 
k = 0.1;
nx = x_bound/dx; ny = y_bound/dy;
dx2 = dx*dx; dy2 = dy*dy;
dt = (dx2 / (2*k));
t_end = 40 * dt;
nt = t_end / dt;

alphax = k*dt / dx2; alphay = k*dt / dy2;

% create mesh
x = 0:dx:x_bound;
y = 0:dy:y_bound;
[X, Y] = meshgrid(x);
% Set Initial Conditions
f = X .* (1-X.^5) .* Y .* (1-Y);
f(1,:) = 0;
f(:,1) = 0;

%% Get Exact Solution
Ntrunc = 50;
f_exact = get_exact(X, Y, Ntrunc, t_end);
f_exact(1,:) = 0;
f_exact(:,1) = 0;

%% Part A) L_inf Error Between f_init, f_exact @ 0

A_exact = get_exact(X, Y, Ntrunc, 0);
L_inf = max(max(abs(f - A_exact)))

%% Solve the Problem Using ADI

% Make The A - matrices
ex = ones(nx-2, 1);
ey = ones(ny-2,1);
A_x = spdiags([-alphax*ex 2*(1+alphax)*ex -alphax*ex], -1:1, nx-2, nx-2);
A_y = spdiags([-alphay*ey 2*(1+alphay)*ey -alphay*ey], -1:1, ny-2, ny-2);

%A_x = full(A_x);
%A_y = full(A_y);

for i=0:nt
    step = mod(i, 2);
    
    if step == 1
        b = (2 + alphay*dy2)*f(2:nx-1,2:ny-1);
        f_inter = A_x \ b;
        b = (2 + alphax*dx2)*f_inter;
        f_save = A_y \ b;
    else
        b = (2 + alphax*dx2)*f(2:nx-1,2:ny-1);
        f_inter = A_y \ b;
        b = (2 + alphay*dy2)*f_inter;
        f_save = A_x \ b;
    end
    f(2:nx-1,2:ny-1) = f_save;
end

figure(1)
surf(X,Y,f)

figure(2)
surf(X,Y,abs(f-f_exact) );

end 



function exactf = get_exact(X, Y, trunc, t_end)
global k
exactf = zeros(length(X), length(Y));
for n=1:trunc
    for m=1:trunc
    exactf=exactf-120*(-n^4*pi^4*(-1)^n + 12*n^2*pi^2*(-1)^n + 24 + 24*(-1)^(1+n))...
        * (-2 + 2*(-1)^m) / (n^7*pi^10*m^3) ...
        * sin(n*pi*X) .*sin(m*pi*Y)*exp(-(n^2+m^2)*pi^2*k*t_end);
    end
end
end
