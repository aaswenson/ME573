%% Alex Swenson ME573 HW10, Problem 1
function problem1_FTCS()
clear; clc;

%% Problem Parameters

dx = 0.1; dy=dx;
dt=0.001;
x = -1:dx:1;
y = 0:dy:2;

nx = length(x); ny = length(y);

% Numerical Parameters
k = 0.5;
nu=0.2;
alphax = nu*dt / dx^2; alphay = nu*dt / dy^2;

% Meshgrid
[X, Y] = meshgrid(x,y);
X=X'; Y=Y';

% % Exact Solution
[u_exact, v_exact] = get_exact(X,Y,k,nu);

%% Initial & Boundary Conditions
fu = X.*Y*0;
fv = X.*Y*0;

fu(:,1) = u_exact(:,1); fu(:,nx) = u_exact(:,nx);
fu(1,:) = u_exact(1,:); fu(ny,:) = u_exact(ny,:);

fv(:,1) = v_exact(:,1); fv(:,nx) = v_exact(:,nx);
fv(1,:) = v_exact(1,:); fv(ny,:) = v_exact(ny,:);

%%
fu_ = fu;
fu_save = fu;

fv_ = fv;
fv_save = fv;

for t=0:dt:15
        Co_x = fu * dt / (dx^2);
        Co_y = fu * dt / (dy^2);
        for j=2:nx-1
            % calculate u direction
            fu_(2:ny-1,j) = fu(2:ny-1,j)-(Co_x(2:ny-1,j)./2).*(fu(3:ny,j)-fu(1:ny-2,j))+...
                alphax*(fu(3:ny,j) - 2*fu(2:ny-1,j) + fu(1:ny-2,j));
            % calculate v direction
            fv_(2:ny-1,j) = fv(2:ny-1,j)-(Co_x(2:ny-1,j)./2).*(fv(3:ny,j)-fv(1:ny-2,j))+...
                alphax*(fv(3:ny,j) - 2*fv(2:ny-1,j) + fv(1:ny-2,j));
        end
        for i=2:ny-1
            % calculate u direction
            fu_save(i,2:nx-1) = fu_(i,2:nx-1)-(Co_y(i,2:nx-1)./2).*(fu_(i,3:nx)-fu_(i,1:nx-2))+...
                alphay*(fu_(i,3:nx) - 2*fu_(i,2:nx-1) + fu_(i,1:nx-2));
            % calculate v direction
            fv_save(i,2:nx-1) = fv_(i,2:nx-1)-(Co_y(i,2:nx-1)./2).*(fv_(i,3:nx)-fv_(i,1:nx-2))+...
                alphay*(fv_(i,3:nx) - 2*fv_(i,2:nx-1) + fv_(i,1:nx-2));
        end
    fu = fu_save;
    fv = fv_save;
end
figure(1)
surf(X,Y,abs(fu-u_exact));
figure(2)
surf(X,Y,abs(fv-v_exact));






end

function [u,v] = get_exact(X, Y, k, nu)

    x0 = 1;
    a0=0.001*k*exp((1+x0)*k); a1=a0; a2=0; a3=a2;

    u_num = -2*nu*(a1 + a3*Y + k*( exp(k*(X-x0) ) - exp(-k*(X-x0)) ).*cos(k*Y));
    v_num = -2*nu*(a2 + a3*X - k*( exp(-k*(X-x0) ) + exp(k*(X-x0)) ).*cos(k*Y));
    
    denom = (a0 + a1.*X +a2*Y + a3*X*Y +(exp(-k*(X-x0)) + exp(k*(X-x0))).*cos(k*Y));
    
    u = u_num ./ denom;
    v = v_num ./ denom;
end