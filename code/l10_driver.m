%% L11 Driver
clc; clear all; close all;

x0 = 1;
xd0 = 0;

[T1,X1] = ode45(@(t,y) l10_lin(t,y),[0,10],[x0 xd0]);
[T2,X2] = ode45(@(t,y) l10_nonlin(t,y),[0,10],[x0 xd0]);


plot(T1,X1(:,1),'r-',T2,X2(:,1),'k-');
set(gca,'YDir','reverse');

grid;axis tight;
legend('Linearized','Non Linear');
xlabel('Time [s]', 'fontsize',14); ylabel('Displacement [m]', 'fontsize',14); 
set(gca,'fontsize',14);