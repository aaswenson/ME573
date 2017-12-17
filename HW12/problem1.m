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
nslip = 10; sliptol = 0.00001; vlid = 1.0; 

alphax = nu * dt / dx^2; alphay = nu * dt / dy^2;
u(Nx,Ny) = 0.0; v(Nx,Ny) = 0.0;

% boundary conditions
utop(1:Nx)=vlid; ubot(1:Nx)=0.0; vleft(1:Ny)=0.0; vright(1:Ny)=0.0;

u(:,Ny) = vlid;

for i=1:tfinal / dt;
    for i=1:nslip;
        [u_1, v_1] = advection_diffusion(dt, dx, dy, u, v, nu,...
                                         ubot', utop', vleft, vright);
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
figure(1)
surf(X,Y,u)
figure(2)
surf(X,Y,v)
end

function [u_2, v_2] = project_vel(Nx, Ny, u, v, dx, dy, dt, rho, phi,...
                                  ubot, utop, vleft, vright)
    u_2 = u; v_2 = v;
    % compute derivative from converged phi (provided)
    dpx = (phi(3:Nx,2:Ny-1) - phi(1:Ny-2,2:Ny-1)) ./ (2*dx);
    dpy = (phi(2:Nx-1,3:Ny) - phi(2:Nx-1,1:Ny-2)) ./ (2*dy);
    % Lower Boundary Derivatives
    Li = (phi(3:Nx-2,1) - phi(2:Nx-3,1)) ./ (2*dx);
    L0 = 0;
    L2 = (-3*phi(2,1) + 4*phi(3,1) - phi(4,1)) ./ (2*dx);
    LN_1 = (3*phi(Nx-1,1) - 4*phi(Nx-2,1) + phi(Nx-3,1)) ./ (2*dx);
    % Upper Boundary Derivatives
    Ui = (phi(3:Nx-2,Ny) - phi(2:Nx-3,Ny)) ./ (2*dx);
    U0 = 0;
    U2 = (-3*phi(2,Ny) + 4*phi(3,Ny) - phi(4,Ny)) ./ (2*dx);
    UN_1 = (3*phi(Nx-1,Ny) - 4*phi(Nx-2,Ny) + phi(Nx-3,Ny)) ./ (2*dx);
    % LHS Boundary Derivatives
    LHSi = (phi(1,3:Ny-2) - phi(1,2:Ny-3)) ./ (2*dy);
    LHS0 = 0;
    LHS2 = (-3*phi(1,2) + 4*phi(1,3) - phi(1,4)) ./ (2*dy);
    LHSN_1 = (3*phi(1,Ny-1) - 4*phi(1,Ny-2) + phi(1,Ny-3)) ./ (2*dy);
    % RHS Boundary Derivatives
    RHSi = (phi(Nx,3:Ny-2) - phi(Nx,2:Ny-3)) ./ (2*dy);
    RHS0 = 0;
    RHS2 = (-3*phi(Nx,2) + 4*phi(Nx,3) - phi(Nx,4)) ./ (2*dy);
    RHSN_1 = (3*phi(Nx,Ny-1) - 4*phi(Nx,Ny-2) + phi(Nx,Ny-3)) ./ (2*dy);

    % project internal velocities
    u_2(2:Nx-1,2:Ny-1) = u(2:Nx-1,2:Ny-1) - (dt/rho)*dpx;
    v_2(2:Nx-1,2:Ny-1) = v(2:Nx-1,2:Ny-1) - (dt/rho)*dpy;
    
    % project bottom nodes
    u_2(3:Nx-2,1) = ubot(3:Nx-2) - (dt/rho).*Li;
    u_2(2,1) = ubot(2) - (dt/rho).*L2;
    u_2(Nx-1,1) = ubot(Nx-1) - (dt/rho).*LN_1;
    v_2(2:Nx-1,1) = 0;
    % project top nodes
    u_2(3:Nx-2,Ny) = utop(3:Nx-2) - (dt/rho).*Ui;
    u_2(2,Ny) = utop(2) - (dt/rho).*U2;
    u_2(Nx-1,1) = ubot(Nx-1) - (dt/rho).*UN_1;
    v_2(2:Nx-1,Ny) = 0;
    % project left nodes
    u_2(1,2:Ny-1) = -(dt/rho)*LHS0;
    v_2(1,3:Ny-2) = vleft(3:Ny-2) - (dt/rho).*LHSi;
    v_2(1,2) = vleft(2) - (dt/rho).*LHS2;
    v_2(1,Ny-1) = vleft(Ny-1) - (dt/rho).*LHSN_1;
    % project right nodes
    u_2(Nx,2:Ny-1) = -(dt/rho)*RHS0;
    v_2(Nx,3:Ny-2) = vright(3:Ny-2) - (dt/rho)*RHSi;
    v_2(Nx,2) = vleft(2) - (dt/rho).*RHS2;    
    v_2(Nx,Ny-1) = vleft(Ny-1) - (dt/rho).*RHSN_1;
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










