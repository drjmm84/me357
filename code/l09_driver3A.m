%% L09 Driver 3A
clc; close all; clear;

[T,X] = ode45(@(t,x) l09_ode3A(t,x),[0:0.01: 60],[0 0 0 0 0 0 0]);

C = [
    1 0 0 0 0 0 0
    0 0 1 0 0 0 0
    0 0 0 0 0 1 0
    ];

Y = (C*X')';    % need transpose from how Matlab saves X

%% Plot 
clf;
plot(T,Y);
legend('x','y','z');
xlabel('Time [s]','fontsize',14);
ylabel('Displacement [m]','fontsize',14);
set(gca,'fontsize',14);
grid on;
axis tight;
