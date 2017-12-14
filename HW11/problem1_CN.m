%% Alex Swenson ME573 HW10, Problem 1
function [X,Y, u_exact, v_exact, fuCN, fvCN] = problem1_CN()
global alphax alphay
%% Problem Parameters

dx = 0.1; dy=dx;
dt=0.1;
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

fuCN = fu;
fvCN = fv;
fu_CN = fu;
fu_saveCN = fu;
fv_CN = fv;    
fv_saveCN = fv;
%% March in Time
for t=0:dt:15
        for j=2:nx-1
            % ----------------- Crank Nicholson ------------------- %
            % calculate u direction
            Aux = update_Ax(fuCN, nx, dt, dx, j);
            bux = update_bx(fuCN, nx, dt, dx, j);
            fu_CN(2:ny-1,j) = thomas(Aux, bux);
            % calculate v direction
            Avx = update_Ax(fvCN, nx, dt, dx, j);
            bvx = update_bx(fvCN, nx, dt, dx, j);
            fv_CN(2:ny-1,j) = thomas(Avx, bvx);
        end
        for i=2:ny-1
            % ----------------- Crank Nicholson ------------------- %
            % calculate u direction
            Auy = update_Ay(fu_CN, ny, dt, dy, i);
            buy= update_by(fu_CN, ny, dt, dy, i);
            fu_saveCN(i,2:nx-1) = thomas(Auy, buy);
            % calculate v direction
            Avy = update_Ay(fv_CN, ny, dt, dy, i);
            bvy= update_by(fv_CN, ny, dt, dy, i);
            fv_saveCN(i,2:nx-1) = thomas(Avy, bvy);
        end
    fuCN = fu_saveCN;
    fvCN = fv_saveCN;
end

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

function A = update_Ax(f, N, dt, dx, j)
    global alphax
    f = f(:,j)';
    e = ones(N-1, 1);
    % set ghost point to make matrix
    f(N+1) = 1;
    A = spdiags([-((f(3:N+1)'*dt/dx) +2*alphax).*e, 4*(1+alphax)*e, (f(1:N-1)'*(dt/dx) - 2*alphax).*e], -1:1, N-2, N-2);
end

function A = update_Ay(f, N, dt, dx, i)
    global alphay
    f = f(i,:);
    e = ones(N-1, 1);
    % set ghost point to make matrix
    f(N+1) = 1;
    A = spdiags([-((f(3:N+1)'*dt/dx) +2*alphay).*e, 4*(1+alphay)*e, (f(1:N-1)'*(dt/dx) - 2*alphay).*e], -1:1, N-2, N-2);
end

function b = update_bx(f, N, dt, dx, j)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    global alphax
    f = f(:,j)';
    c = dt/dx;
    b = (c*f(2:N-1)+2*alphax).*f(1:N-2) + 4*(1-alphax)*f(2:N-1) + (2*alphax-c*f(2:N-1)).*f(3:N);
    b(1) = b(1) + (c*f(2)+2*alphax)*f(1);
    b(N-2) = b(N-2) + (2*alphax - f(N-2)*c)*f(N);
end

function b = update_by(f, N, dt, dy, i)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    global alphay
    f = f(i,:);
    c = dt/dy;
    b = (c*f(1:N-2)+2*alphay).*f(1:N-2) + 4*(1-alphay)*f(2:N-1) + (2*alphay-c*f(2:N-1)).*f(3:N);
    b(1) = b(1) + (c*f(2)+2*alphay)*f(1);
    b(N-2) = b(N-2) + (2*alphay - f(N-2)*c)*f(N);
end



