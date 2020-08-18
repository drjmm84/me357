%% Example 1
clc; clear; close all;

t = linspace(0,5,1e3);      % create time vector
x = exp(-2*t)-exp(-3*t);    % create displacement
% plot
plot(t,x,'k-');
grid on;
xlabel('Time [s]','fontsize',14);
ylabel('Displacement [m]','fontsize',14);
set(gca,'fontsize',14);

%% Example 2
clc; clear; close all;

t = linspace(0,10,1e3);         % create time vector
x = exp(-t)/5 - exp(-6*t)/5;    % create displacement vector
% plot
plot(t,x,'r-');
grid on;
xlabel('Time [s]','fontsize',14);
ylabel('Displacement [m]','fontsize',14);
set(gca,'fontsize',14);

%% Together
clc; clear; close all;

t = linspace(0,6,2e3);

x1 = exp(-2*t)-exp(-3*t);       % system 1
x2 = exp(-t)/5 - exp(-6*t)/5;   % system 2

plot(t,x1,'k-',t,x2,'r-');
grid on;
xlabel('Time [s]','fontsize',14);
ylabel('Displacement [m]','fontsize',14);
legend('c=5 N-s/m','c=7 N-s/m','location','northeast');
set(gca,'fontsize',14);