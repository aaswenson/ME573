%% Alex Swenson ME573 HW10, Problem 1
function [u_2, v_2] = advection_diffusion(dt, dx, dy, u, v, nu,...
                                           ubot, utop, vleft, vright)
global alphax alphay
%% Problem Parameters

nx = length(u(:,1));
ny = length(v(1,:));
alphax = nu*dt / dx^2; alphay = nu*dt / dy^2;

u_1 = u; u_2 = u;
v_1 = v; v_2 = v;

%% March in Time
    for j=2:nx-1
        % ----------------- Crank Nicholson ------------------- %
        % calculate u direction
        Aux = update_Ax(u, nx, dt, dx, j);
        bux = update_bx(u, u, nx, dt, dx, j, 0, 0);
        u_1(2:ny-1,j) = thomas(Aux, bux);
        % calculate v direction
        Avx = update_Ax(v, nx, dt, dx, j);
        bvx = update_bx(u, v, nx, dt, dx, j, vleft(j), vright(j));
        v_1(2:ny-1,j) = thomas(Avx, bvx);
    end
    for i=2:ny-1
        % ----------------- Crank Nicholson ------------------- %
        % calculate u direction
        Auy = update_Ay(u_1, ny, dt, dy, i);
        buy= update_by(v, u_1, ny, dt, dy, i, ubot(i), utop(i));
        u_2(i,2:nx-1) = thomas(Auy, buy);
        % calculate v direction
        Avy = update_Ay(v_1, ny, dt, dy, i);
        bvy= update_by(v, v_1, ny, dt, dy, i, 0, 0);
        v_2(i,2:nx-1) = thomas(Avy, bvy);
    end


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

function b = update_bx(u, f, N, dt, dx, j, bc1, bc2)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    global alphax
    f = f(:,j)';
    u = u(:,j)';
    c = dt/dx;
    b = (c*u(2:N-1)+2*alphax).*f(1:N-2) + 4*(1-alphax)*f(2:N-1) + (2*alphax-c*u(2:N-1)).*f(3:N);
    b(1) = b(1) + (c*u(2)+2*alphax)*bc1;
%     b(N-2) = b(N-2) + (2*alphax - u(N-2)*c)*bc2;
end

function b = update_by(v, f, N, dt, dy, i, bc1, bc2)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    global alphay
    f = f(i,:);
    v = v(i,:);
    c = dt/dy;
    b = (c*v(1:N-2)+2*alphay).*f(1:N-2) + 4*(1-alphay)*f(2:N-1) + (2*alphay-c*v(2:N-1)).*f(3:N);
    b(1) = b(1) + (c*v(2)+2*alphay)*bc1;
%     b(N-2) = b(N-2) + (2*alphay - v(N-2)*c)*bc2;
end



