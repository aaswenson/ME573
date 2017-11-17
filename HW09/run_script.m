clear; clc;

disp('Running FTCS and FTBS algorithms...')
[~, ~, f_ftbs, f_ftcs] = problem2(0.01, 0.05, 0.5);
disp('Running FTCS and FTBS in unstable configuration')
[x_unstable, exact_unstable, ftbs_unstable, ftcs_unstable] = problem2(0.02, 0.025, 0.5);
disp('Running CN algorithm...')
[x, f_analytic, f_CN] = problem3(0.01, 0.05, 0.5);


%% Calculate Error
CN_spatial_error = abs(f_CN - f_analytic);
FTBS_spatial_error = abs(f_ftbs - f_analytic);


%% Plotting
start_plot = 1;
i = 0;
% Plot Results 
figure(start_plot + i)
plot(x, f_analytic, 'r', x, f_ftcs, 'g', x, f_ftbs, 'b', x, f_CN, 'y')
title('Comparison of Linear Advection Schemes')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('Analytic', 'FTCS', 'FTBS', 'CN')
saveas(gcf,'./writeup/ftcs_ftbs.png')
i = i + 1;

% Plot Unstable Results
figure(start_plot + i)
plot(x_unstable, exact_unstable, 'r', x_unstable, ftbs_unstable, 'b', x_unstable, ftcs_unstable, 'g')
title('Comparison of Linear Advection Schemes in an Unstable State')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('Analytic', 'FTBS', 'FTCS')
saveas(gcf,'./writeup/ftcs_ftbs_unstable.png')
i = i + 1;

% Plot results without FTCS
figure(start_plot + i)
plot(x, f_analytic, 'r', x, f_ftbs, 'b', x, f_CN, 'g')
title('Comparison of Linear Advection Schemes')
xlabel('X [-]')
ylabel('f(X) [-]')
legend('Analytic', 'FTBS', 'CN')
saveas(gcf,'./writeup/ftbs_CN.png')
i = i + 1;

% Plot Spatial Error For CN and FTBS
figure(start_plot + i)
plot(x, CN_spatial_error, 'r', x, FTBS_spatial_error, 'g')
title('Spatial Error')
xlabel('X [-]')
ylabel('Error [-]')
legend('CN', 'FTBS')
saveas(gcf,'./writeup/ftbs_CN_error.png')
i = i + 1;