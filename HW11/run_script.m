%% HW11 Run Script
clc; clear; clf;
tic;
[ftcsx, ftcsy, ftcsu, ftcsv, ftcsuexact, ftcsvexact] = problem1_FTCS();
t_ftcs = toc;
tic;
[CNx, CNy, CNuexact, CNvexact, CNu, CNv] = problem1_CN();
t_CN = toc;

disp(['The FTCS method took ', num2str(t_ftcs), ' seconds.']);
disp(['The CN method took ', num2str(t_CN), ' seconds.']);

error_CN_u = abs(CNuexact - CNu);
error_CN_v = abs(CNvexact - CNv);
error_ftcs_u = abs(ftcsuexact -ftcsu);
error_ftcs_v = abs(ftcsvexact -ftcsv);

%% Plotting
start_plot = 1;
i = 0;
% Plot Results 
figure(start_plot + i)
surf(ftcsx, ftcsy, error_ftcs_u)
title('FTCS U Error')
xlabel('X [-]')
ylabel('f(X) [-]')
saveas(gcf,'./writeup/ftcsu.png')
i = i + 1;

% Plot Results 
figure(start_plot + i)
surf(ftcsx, ftcsy, error_ftcs_v)
title('FTCS V Error')
xlabel('X [-]')
ylabel('f(X) [-]')
saveas(gcf,'./writeup/ftcsv.png')
i = i + 1;

% Plot Results 
figure(start_plot + i)
surf(CNx, CNy, error_CN_u)
title('Crank Nicolson U Error')
xlabel('X [-]')
ylabel('f(X) [-]')
saveas(gcf,'./writeup/CNu.png')
i = i + 1;

% Plot Results 
figure(start_plot + i)
surf(CNx, CNy, error_CN_v)
title('Crank Nicolson V Error')
xlabel('X [-]')
ylabel('f(X) [-]')
saveas(gcf,'./writeup/CNv.png')
i = i + 1;