%% L09 1A
clc; close all; clear;

[T,X] = ode45(@(t,x) l09_ode1A(t,x),[0:0.01:10],[-1 0]);    % solve SS eqs 

C = [1 0];
D = [0];        %#ok<*NBRAK>

Y = (C*X')';    % need transpose from how Matlab saves X

%% Plot
plot(T,Y,'k-');
% legend('x');
xlabel('Time [s]','fontsize',14);
ylabel('Displacement [m]','fontsize',14);
set(gca,'fontsize',14);
grid on;
