%% HW10 Run Script
clc; clear; clf;
[p1x, p1ftbs, p1ftcs, p1CN, p1exact] = problem1();
[p2x, p2ftcs, p2CN, p2exact] = problem2();


%% Plotting
start_plot = 1;
i = 0;
% Plot Results 
figure(start_plot + i)
plot(p1x, p1ftbs, 'r', p1x, p1ftcs, 'g', p1x, p1CN, 'b', p1x, p1exact, '-o')
title('Comparison of Linear Advection Schemes')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('FTBS', 'FTCS', 'CN', 'Exact')
saveas(gcf,'./writeup/problem1.png')
i = i + 1;

figure(start_plot + i)
plot(p2x, p2ftcs, 'o', p2x, p2exact)
title('Burgers Equation')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('FTCS', 'Exact')
saveas(gcf,'./writeup/problem2_ftcs.png')
i = i + 1;

figure(start_plot + i)
plot(p2x, p2CN, 'o', p2x, p2exact)
title('Burgers Equation')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('Crank-Nicolson', 'Exact')
saveas(gcf,'./writeup/problem2_CN.png')
i = i + 1;
