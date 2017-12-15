%Validation data from
%U. GHIA, K. N. GHIA, AND C. T. SHIN (1982) "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a
%Multigrid Method,"  JOURNAL OF COMPUTATIONAL PHYSICS 48, 387-411
%
clear; clc ; close all;
load ('N50_check')   % load data from computation
midpoint=26;   % Midpoint where the data resides
%Table 8 Tabulated u-velocity profiles along a vertical
%line passing through the geometric center of the cavity at variable Re.
y_p=[1.0000 0.9766  0.9688 0.9609  0.9531 0.8516  0.7344 0.6172 0.5000 0.4531 0.2813 0.1719  0.1016 0.0703 0.0625 0.0547 0.0000];%y coordinate
u_re100=[1.0000 0.8412 0.7887 0.7372 0.68717 0.2315 0.0033  -0.1364  -0.2058  -0.2109  -0.1566  -0.1015  -0.0643  -0.04775  -0.0419  -0.0371 0.0000];% Re=100
%Table 9 Tabulated v-velocity profiles along a horizontal line passing through the geometric center of the cavity at various Reynolds numbers.

x_p=[1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5000 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0.0000];
v_re100=[0.0000 -0.05906  -0.0739 -0.0886 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.1009 0.0923 0.0000];
% Post-processing
%
figure('units','normalized','position',[0.0 0.5 .4 .4])
plot(y_p,u_re100,'-ro',Y(1,:),u(midpoint,:),'-+')
grid on
set(gca,'fontsize',18)
legend('Re=100 data','Re=100 computation')
title('u vs Y comparison (at X=0.5) ')
xlabel('y')
ylabel('u')
%
figure('units','normalized','position',[0.5 0.5 .4 .4])
plot(x_p,v_re100,'-ro',X(:,1),v(:,midpoint),'-+')
grid on
set(gca,'fontsize',18)
legend('Re=100 data','Re=100 computation')
title('v vs X comparison (at Y=0.5) ')
xlabel('x')
ylabel('v')
