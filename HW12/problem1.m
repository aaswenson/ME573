%% HW12 Alexander Swenson
function problem1()
%% Set Parameters

nu = 0.01; rho = 1.0; dt = 0.1; tfinal = 15;
nslip = 50; sliptol = 0.00001; dx = 0.02; dy=dx; Lx=1.0; Ly=Lx;
Nx = Lx / dx;   Ny = Ly / dy;
global alphax alphay
alphax = nu * dt / dx^2; alphay = nu * dt / dy^2;

u(Nx,Ny) = 0.0; v(Nx,Ny) = 0.0;

end

function A = update_Ax(f, N, dt, dx, j)
    global alphax
    f = f(:,j)';
    e = ones(N-1, 1);
    % set ghost point to make matrix
    f(N+1) = 1;
    A = spdiags([-((f(3:N+1)'*dt/dx) +2*alphax).*e, 4*(1+alphax)*e, (f(1:N-1)'*(dt/dx) - 2*alphax).*e], -1:1, N-2, N-2);
end