%% Run Problem 1
clear; clc; clear; clf;
dx_vals = [0.01 0.05 0.1];

for i=1:length(dx_vals)
problem1(dx_vals(i),i);    
end
