%% Lecture 12 Example 1
clear; close all; clc; 

b = 50;                         % step input amplitude [VDC]
% Create TFs 
G1 = tf([1],[2 1]);             %#ok<*NBRAK>
G2 = tf([1],[3 1]);

% System 1
figure(1);
step(b*G1,0:0.1:30);            % apply step input of 50 VDC
grid on;
ylabel('$v_1$ [VDC]','fontsize',14,'interpreter','LaTex');
set(gca,'fontsize',14);

% System 2
figure(2);
step(b*G2,0:0.1:30);            % apply step input of 50 VDC
grid on;
ylabel('$v_1$ [VDC]','fontsize',14,'interpreter','LaTex');
set(gca,'fontsize',14);

% Both systems 
figure(3);
step(b*G1,0:0.1:30); hold on;   % plot two systems on one plot 
step(b*G2,0:0.1:30);
grid on;
ylabel('$v_1$ [VDC]','fontsize',14,'interpreter','LaTex');
legend('System 1','System 2','location','southeast');
set(gca,'fontsize',14);

%% Lecture 12 Example 2
clear; close all; clc; 

G = tf([1.5 1],[3.6 1]);        % create TF
b = 5;                          % step amplitude [VDC]
[y,t]=step(b*G,0:0.01:25);
plot(t,y,'b-')
axis([0 25 0 5.2]);
grid on;
ylabel('$v_o$ [VDC]','fontsize',14,'interpreter','LaTex');
xlabel('Time [s]','fontsize',14);

set(gca,'fontsize',14);
