%% L10 Driver: Air Resistance 
clc; clearvars; close all;

x0 = -50e-2;
xd0 = 0;
tf = 20;
[T1,X1] = ode45(@(t,y) l10A_lin(t,y),[0,tf],[x0 xd0]);
[T2,X2] = ode45(@(t,y) l10A_nonlin(t,y),[0,tf],[x0 xd0]);


plot(T1,1e2*X1(:,1),'r-',T2,1e2*X2(:,1),'k-');
set(gca,'YDir','reverse');

grid;   axis tight;
legend('Linearized','Non Linear');
xlabel('Time [s]', 'fontsize',14); ylabel('Displacement [cm]', 'fontsize',14); 
set(gca,'fontsize',14);