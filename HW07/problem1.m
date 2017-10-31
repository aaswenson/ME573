%% Alex Swenson HW7 (Computer Project 5c)
function problem1()

%% Set Domain

% some important constants 
x_bound = 1; y_bound = 1.;
dx = 0.05; dy = dx;
global k 
k = 0.1;
nx = x_bound/dx; ny = y_bound/dy;
dx2 = dx*dx; dy2 = dy*dy;
dt = (dx2 / (2*k));
t_end = 40 * dt;
nt = t_end / dt;

alphax = k*dt / dx2; alphay = k*dt / dy2;

x = 0:dx:x_bound;
y = 0:dy:y_bound;
% create mesh
[X, Y] = meshgrid(x);
% Set Initial Conditions
f = X .* (1-X.^5) .* Y .* (1-Y);

%% Get Exact Solution
Ntrunc = 50;
f_exact = get_exact(X, Y, Ntrunc, t_end);

%% Part A) L_inf Error Between f_init, f_exact @ 0

A_exact = get_exact(X, Y, Ntrunc, 0);

figure(1)
surf(X, Y, abs(f - A_exact));

% Make The A - matrices

e = ones(nx-2, 1);
A_x = spdiags([alphax*e 2*(1+alphax)*e -alphax*e], -1:1, nx-2, nx-2);
A_y = spdiags([alphay*e 2*(1+alphay)*e -alphay*e], -1:1, ny-2, ny-2);

A_x = full(A_x);
A_y = full(A_y);

f_save = f;

for i=0:nt
    step = mod(i, 2);
    
    if step == 0
        for j=2:ny-2
       % use the first routine (x then y)
            f_intermediate =  A_x \ f(j,1:ny-2)';
            f_save(j,1:ny-2) = (A_y \ f_intermediate)';
        end
    else
        for j=2:nx-2
       % use the second solving routine (y then x)
            f_intermediate = A_y \ f(j,1:nx-2)';
            f_save(j,1:nx-2) = (A_x \ f_intermediate)';
        end
    f(1:ny-2,1:ny-2) = f_save
    end
end

figure(2)
surf(X,Y,f)

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
