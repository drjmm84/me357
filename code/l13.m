%% Lecture 13

%% Example 1

clear; close all; clc; 
figure(1);
G = tf([1],[5 10 500]);         %#ok<*NBRAK>
[y, t]=step(10*G,0:0.001:5);     % simulate step function input of amplitude 10 
% Plot 
plot(t,1000*y,'b-');            % plot in mm
grid on;
ylabel('y [mm]','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

%% Example 2
clear; close all; clc; 

z = 0.15;                       % damping ratio
wn = 25;                        % natural freq. [rad/s]
yss = 24;                       % steady-state value 
wd = wn*sqrt(1-z^2);            % damped freq. [rad/s]
T = 2*pi/wd;                    % period [s]
disp(T*[0.5 1 1.5 2 2.5 3]);    % location of peaks/valleys
y0 = yss*(1+exp(-z*wn*1*T/2));  % max value / first peak
y05 = yss*(1-exp(-z*wn*1*T));   % first valley
y1 = yss*(1+exp(-z*wn*1.5*T));  % second peak
y15 = yss*(1-exp(-z*wn*2*T));   % second valley

figure(1);
G = tf([8*wn^2],[1 2*z*wn wn^2]);
[y, t]=step(3*G,0:0.0001:1.4);
plot(t, y, 'b-'); hold on;
grid on;
plot(t,8*3*(1+exp(-z*wn*t)),'k--'); % bounding envelopes 
plot(t,8*3*(1-exp(-z*wn*t)),'k--');
axis tight;
ylabel('y [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

%% Example 3: SS model

clear; close all; clc; 
figure(1);

A = [0 1 
    -10 -2.5];

L = eig(A);
disp(L);            % eigenvals/poles of TF
if ~isreal(L)       % check if complex
    wn = abs(L);    % natural freq. [rad/s]
    disp(wn);
end

% convert to TF
[b,a] = ss2tf(A,[0;0.5],[1 0],0);
disp(roots(a));     % poles of TF
G = tf(b,a); 
[y,t]=step(10*G,0:0.01:5);
plot(t,y,'b-'); hold on;
grid on;
% axis tight;
ylabel('y [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

%% Extra Pole
clc; close all; clear;
roots([1 4 24.542])                             % dominant poles
G1 = tf(24.542,[1 4 24.542]);                   % base TF
G2 = tf(24.542,[0.1 1.4 6.4542 24.542]);        % extra pole at -10
G3 = tf(24.542,[0.3333 2.3333 12.1807 24.542]); % extra pole at -3 
[Y1,T] = step(G1,0:0.01:3.5);
[Y2,~] = step(G2,T);
[Y3,~] = step(G3,T);
plot(T,Y1,'r-',T,Y2,'b-',T,Y3,'g-');
legend('G_1','G_2','G_3','location','southeast');
grid on; set(gca,'fontsize',14);
ylabel('x [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
axis tight;

%% Extra Zero
clc; close all; clear;

G1 = tf(24.542,[1 4 24.542]);               % base TF
G2 = tf([2.4542 24.542],[1 4 24.542]);      % extra zero at -10 
G3 = tf([8.18067 24.542],[1 4 24.542]);     % extra zero at -3
[Y1,T] = step(G1,0:0.01:3);
[Y2,~] = step(G2,T);
[Y3,~] = step(G3,T);
plot(T,Y1,'r-',T,Y2,'b-',T,Y3,'g-');
legend('G_1','G_2','G_3','location','southeast');
grid on; set(gca,'fontsize',14);
ylabel('x [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
axis tight;

%% Extra RHP Zero
clc; close all; clear;

G1 = tf(24.542,[1 4 24.542]);               % base TF
G2 = tf([-12.271 24.542],[1 4 24.542]);     % extra zero at +2
G3 = tf([-2.72689 24.542],[1 4 24.542]);    % extra zero at +9
[Y1,T] = step(G1,0:0.01:3.5);
[Y2,~] = step(G2,T);
[Y3,~] = step(G3,T);
plot(T,Y1,'r-',T,Y2,'b-',T,Y3,'g-');
legend('G_1','G_2','G_3','location','southeast');
grid on; set(gca,'fontsize',14);
ylabel('x [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
axis tight;

%% Extra Zero Close to Extra Pole 
clc; close all; clear;

G1 = tf(24.542,[1 4 24.542]);                       % base TF
% G2 = tf([26.1781 78.5344],[1 7.2 37.342 78.5344]);
% G2 = tf([8.181 24.542],[0.3125 2.25 11.67 24.54]);
G2 = tf([24.542],[0.08 1.32 5.963 24.542]);         % p = -12.5
G3 = tf([2.028 24.542],[0.08 1.32 5.963 24.542]);   % p = -12.5, z = -12.1
[Y1,T] = step(G1,0:0.01:3);
[Y2,~] = step(G2,T);
[Y3,~] = step(G3,T);
plot(T,Y1,'r-',T,Y2,'b-',T,Y3,'g-');
legend('G_1','G_2','G_3','location','southeast');
grid on; set(gca,'fontsize',14);
ylabel('x [m]','fontsize',14);
xlabel('Time [s]','fontsize',14);
axis tight;
