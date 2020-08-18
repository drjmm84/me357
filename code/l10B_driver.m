%% L10 Driver
clc; clearvars; close all;

global F;
F = 2;

x0 = 0;
xd0 = 0;
tspn = 200;

[T1,X1] = ode45(@(t,y) l10B_lin(t,y),[0,tspn],[x0 xd0]);
[T2,X2] = ode45(@(t,y) l10B_nonlin(t,y),[0,tspn],[x0 xd0]);


plot(T1,X1(:,1),'r-',T2,X2(:,1),'k-');
set(gca,'YDir','reverse');

grid;axis tight;
legend('Linearized','Non Linear');
xlabel('Time [s]', 'fontsize',14); ylabel('Displacement [m]', 'fontsize',14); 
set(gca,'fontsize',14);
%%
w=1;
f = 2*pi/w;
pc1 = interp1(T2,X2(:,1),[1:30]*f,'pchip');
pc2 = interp1(T2,X2(:,2),[1:30]*f,'pchip');