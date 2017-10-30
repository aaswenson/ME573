%% Alex Swenson HW7 (Computer Project 5c)
function problem1()
clc; clear;
%% Set Domain

% some important constants 
x_bound = 1; y_bound = 1.;
dx = 0.05; dy = dx;
global k 
k = 0.1;
nx = x_bound/dx; ny = y_bound/dy;
dx2 = dx*dx; dy = dy*dy;
dt = (dx2 / k) / 4.0;
t_end = 40 * dt;

x = 0:dx:x_bound;
y = 0:dy:y_bound;
% create mesh
[X, Y] = meshgrid(x);
% Set Initial Conditions
f_init = X .* (1-X.^5) .* Y .* (1-Y);

%% Get Exact Solution
Ntrunc = 50;
f_exact = get_exact(X, Y, Ntrunc, t_end);

%% Part A) L_inf Error Between f_init, f_exact @ 0

A_exact = get_exact(X, Y, Ntrunc, 0);

figure(1)
surf(X, Y, abs(f_init - A_exact));

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
