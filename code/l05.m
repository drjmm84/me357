%% Example 1A
clc; clear; close all;
m = 0.1;
k = 6;
b = 0.4;

G = tf([1],[m b k]); %#ok<*NBRAK>
% Plot 
impulse(10*G); % step function of 10 units
grid on;
ylabel('Displacement [m]','fontsize',14);
xlabel('Time','fontsize',14);
title('');
set(gca,'fontsize',14);

%% Example 1B
clc; clear; close all;
m = 0.1;
k1 = 6;
k2 = 4;
b2 = 0.4;

G = tf([b2 k2],[b2*m k2*m (b2*k1+b2*k2) k1*k2]);

% Step
figure(1);
step(10*G); % step function of 10 N 
grid on;
ylabel('Displacement [m]','fontsize',14);
xlabel('Time','fontsize',14);
title('');
set(gca,'fontsize',14);

% Impulse
figure(2);
impulse(10*G);  % impusle of 10 units
grid on;
ylabel('Displacement [m]','fontsize',14);
xlabel('Time','fontsize',14);
title('');
set(gca,'fontsize',14);
