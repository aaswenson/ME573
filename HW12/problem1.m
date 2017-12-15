%% HW12 Alexander Swenson
function problem1()
clear; clc; clf;
%% Set Parameters
global alphax alphay
S = load('converged_phi.mat');
phi = S.phi;
dx = 0.02; dy=dx; Lx=1.0; Ly=Lx;
Nx = (Lx/dx) + 1;   Ny = (Ly/dy) + 1;
save_idx = [int8(0.14/dx), int8(0.08/dy)];
nu = 0.01; rho = 1.0; dt = 0.1; tfinal = 1; t=0:dt:tfinal;
nslip = 50; sliptol = 0.00001; vlid = 1.0; 

alphax = nu * dt / dx^2; alphay = nu * dt / dy^2;
u(Nx,Ny) = 0.0; v(Nx,Ny) = 0.0;

% boundary conditions
utop(1:Nx)=vlid; ubot(1:Nx)=0.0; vleft(1:Ny)=0.0; vright(1:Ny)=0.0;

u(:,Ny) = vlid;

for i=1:tfinal / dt;
    for i=1:nslip;
        [u_1, v_1] = advection_diffusion(dt, dx, dy, u, v, nu);
        [u_2, v_2] = project_vel(Nx, Ny, u_1, v_1, dx, dy, dt, rho, phi,...
                         ubot', utop', vleft, vright);
        [utop, ubot, vleft, vright, maxslip] =...
            check_vel(u_2, v_2, Nx, Ny, vlid, utop, vleft, ubot, vright);
        if maxslip < sliptol
           disp(['required slip iterations ', num2str(i)])
           break 
        end
    end     
     u = u_2;
     v = v_2;
end

x = 0:dx:1;
y = 0:dy:1;

% Meshgrid
[X, Y] = meshgrid(x,y);
X=X'; Y=Y';

surf(X,Y,u)
end

function [u_2, v_2] = project_vel(Nx, Ny, u, v, dx, dy, dt, rho, phi,...
                                  ubot, utop, vleft, vright)
    % compute derivative from converged phi (provided)
    dpx = (phi(3:Nx,2:Ny-1) - phi(1:Ny-2,2:Ny-1)) ./ (2*dx);
    dpy = (phi(2:Nx-1,3:Ny) - phi(2:Nx-1,1:Ny-2)) ./ (2*dy);
    % Lower Boundary Derivatives
    Ldpdx_int = (phi(3:Nx,1) - phi(1:Nx-2,1)) ./ (2*dx);
    Ldpdy = 0;
    % Upper Boundary Derivatives
    Udpdx_int = (phi(3:Nx,Ny) - phi(1:Nx-2,Ny)) ./ (2*dx);
    Udpdy = 0;
    % LHS Boundary Derivatives
    LHSdpdy_int = (phi(1,3:Ny) - phi(1,1:Ny-2)) ./ (2*dy);
    LHSdpdx = 0;
    % RHS Boundary Derivatives
    RHSdpdy_int = (phi(Nx,3:Ny) - phi(Nx,1:Ny-2)) ./ (2*dy);
    RHSdpdx = 0;

    % project internal velocities
    u_2(2:Nx-1,2:Ny-1) = u(2:Nx-1,2:Ny-1) - (dt/rho)*dpx;
    v_2(2:Nx-1,2:Ny-1) = v(2:Nx-1,2:Ny-1) - (dt/rho)*dpy;
    % project bottom nodes
    u_2(2:Nx-1,1) = ubot(2:Nx-1) - (dt/rho)*Ldpdx_int;
    v_2(2:Nx-1,1) = -(dt/rho) * Ldpdy;
    % project top nodes
    u_2(2:Nx-1,Ny) = utop(2:Nx-1) - (dt/rho)*Udpdx_int;
    v_2(2:Nx-1,Ny) = -(dt/rho) * Udpdy;
    % project left nodes
    u_2(1,2:Ny-1) = -(dt/rho)*LHSdpdx;
    v_2(1,2:Ny-1) = vleft(2:Ny-1) - (dt/rho)*LHSdpdy_int;
    % project right nodes
    u_2(Nx,2:Ny-1) = -(dt/rho)*RHSdpdx;
    v_2(Nx,2:Ny-1) = vright(2:Ny-1) - (dt/rho)*RHSdpdy_int;   
end

function [utop, ubot, vleft, vright, maxslip] = check_vel(u, v, Nx, Ny,...
                                          vlid, utop, vleft, ubot, vright)

    maxslip = 0;
    for i=1:Nx
       % Top
       slip = u(i,Ny) - vlid; utop(i) = utop(i) - slip;
       maxslip = check_slip(maxslip, slip);
       % Bottom
       slip = u(i,1); ubot(i) = ubot(i) - slip;
       maxslip = check_slip(maxslip, slip);
    end
    for j=1:Ny
       % Left
       slip = v(1,j); vleft(j) = vleft(j) - slip; 
       maxslip = check_slip(maxslip, slip);
       % Right
       slip = v(Nx,j); vright(j) = vright(j) - slip;
       maxslip = check_slip(maxslip, slip);
    end
end

function maxslip = check_slip(max, slip)
    maxslip = max;
    if abs(slip) > max
        maxslip = abs(slip);
    end
end










