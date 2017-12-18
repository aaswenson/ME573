%% Alex Swenson ME573 HW10, Problem 1
function [u_2, v_2] = advection_diffusion(dt, dx, dy, u, v, nu,...
                                           ubot, utop, vleft, vright)
%% Problem Parameters

nx = length(u(:,1));
ny = length(v(1,:));
alphax = nu*dt / dx^2; alphay = nu*dt / dy^2;

u_1 = u; u_2 = u;
v_1 = v; v_2 = v;

Cox = u.*(dt/dx); Coy = v.*(dt/dy);
%% March in Time
    for j=2:nx-1
        % ----------------- Crank Nicholson ------------------- %
        % calculate u direction
        Ax = update_A(Cox(:,j)', nx, alphax);
        bux = update_b(Cox(:,j)', u(:,j)', nx, alphax, 0, 0);
        u_1(2:ny-1,j) = thomas(Ax, bux);
        % calculate v direction
        bvx = update_b(Cox(:,j)', v(:,j)', nx, alphax, vleft(j), vright(j));
        v_1(2:ny-1,j) = thomas(Ax, bvx);
    end
    for i=2:ny-1
        % ----------------- Crank Nicholson ------------------- %
        % calculate u direction
        Ay = update_A(Coy(:,j)', ny, alphay);
        buy= update_b(Coy(i,:), u_1(i,:), ny, alphay, ubot(i), utop(i));
        u_2(i,2:nx-1) = thomas(Ay, buy);
        % calculate v direction
        bvy= update_b(Coy(i,:), v_1(i,:), ny, alphay, 0, 0);
        v_2(i,2:nx-1) = thomas(Ay, bvy);
    end


end

function A = update_A(f, N, alpha)
    e = ones(N-1, 1);
    % set ghost point to make matrix
    f(N+1) = 1;
    A = spdiags([-(f(3:N+1)' +2*alpha).*e, 4*(1+alpha)*e, (f(1:N-1)' - 2*alpha).*e], -1:1, N-2, N-2);
end


function b = update_b(Co, f, N, alpha, bc1, bc2)
    % update the RHS vector, include initial boundary conditions for 
    % first and last nodes.
    b = (Co(2:N-1)+2*alpha).*f(1:N-2) + 4*(1-alpha)*f(2:N-1) + (2*alpha-Co(2:N-1)).*f(3:N);
    b(1) = b(1) + (Co(2)+2*alpha)*bc1;
    b(N-2) = b(N-2) + (2*alpha - Co(N-2))*bc2;
end



