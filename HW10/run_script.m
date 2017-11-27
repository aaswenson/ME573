%% HW10 Run Script
clc; clear; clf;
[p1x, p1ftbs, p1ftcs, p1CN, p1exact] = problem1();


%% Plotting
start_plot = 1;
i = 0;
% Plot Results 
figure(start_plot + i)
plot(p1x, p1ftbs, 'r', p1x, p1ftcs, 'g', p1x, p1CN, 'b', p1x, p1exact, '-o')
title('Comparison of Linear Advection Schemes')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('FTBS', 'FTCS', 'CN','Exact')
saveas(gcf,'./writeup/problem1.png')
i = i + 1;