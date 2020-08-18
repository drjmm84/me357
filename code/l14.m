%% Lecture 14

%% Example 0

clc; clear; close all;
wn = 10;
hold on;
for z=0.1:0.1:0.7
    G = tf(1,[1, 2*z*wn, wn^2]);  % create TF
    bode(G,1e-0:0.1:1e2);       % makes Bode Plot with specified freq. resolution
end
legend('\zeta=0.1','\zeta=0.2','\zeta=0.3','\zeta=0.4','\zeta=0.5',...
    '\zeta=0.6','\zeta=0.7','location','northeast');
grid on;
set(gca,'fontsize',13);

%% Example 1

clc; clear; close all;

G = tf(1,[1,10]);       % create TF
bode(G,1e-1:0.1:1e3);   % makes Bode Plot with specified freq. resolution 
grid on;
set(gca,'fontsize',13);

%% Plot Example 1

clc; close all; clear;
A = 10;                         %#ok<*NASGU> % input amplitude 
w = 10;                         % input frequency [rad/s]
sim('l14_1');                   % run Simulink file with this A and w

tim = simout.time(400:end);     % remove transient response 
tim = tim-min(tim);             % rezero time

u = simout.data(400:end,1);     % input data
y = simout.data(400:end,2);     % output data
plot(tim,u,'r',tim,y,'k');
grid on; axis tight;
set(gca,'fontsize',14)
xlabel('Time [s]','fontsize',14);
ylabel('Amplitude','fontsize',14);
legend('Input','Output');

%% Example 2A

clc; clear; close all;

M = 1000;
b = 1000;
k = 10^7;
wn = sqrt(k/M);         % natural freq. [rad/s]
z = b/(2*M*wn);         % damping ratio
wr = wn*sqrt(1-2*z^2);    % resonant freq. [rad/s]

G = tf(1,[M b k]);
w = logspace(1,3,1e4);     % driving freq range

[M, phi] = bode(G,w); grid on;
M = M(:);   phi = phi(:);

subplot(2,1,1);
title('Bode Diagram','fontsize',14);
semilogx(w,20*log10(M),'b-');
grid on;
ylabel('M [dB]','fontsize',14);
set(gca,'fontsize',14);

subplot(2,1,2);
semilogx(w,phi,'b-');
grid on;
ylabel('\phi [\circ]','fontsize',14);
xlabel('\omega [rad/s]','fontsize',14);
set(gca,'fontsize',14);
axis tight;

%% Example 2B: double spring constant 

clc; clear; close all;

M = 1000;
b = 1000;
k = 2*10^7;
wn = sqrt(k/M);         % natural freq. [rad/s]
z = b/(2*M*wn);         % damping ratio
wr = wn*sqrt(1-2*z^2);  % resonant freq. [rad/s]

G = tf(1,[M b k]);
w = logspace(1,3,1e4);     % driving freq range

[M,phi] = bode(G,w); grid on;
M = M(:);   phi = phi(:);

subplot(2,1,1);
title('Bode Diagram','fontsize',14);
semilogx(w,20*log10(M),'b-');
grid on;
ylabel('M [dB]','fontsize',14);
set(gca,'fontsize',14);

subplot(2,1,2);
semilogx(w,phi,'b-');
grid on;
ylabel('\phi [\circ]','fontsize',14);
xlabel('\omega [rad/s]','fontsize',14);
set(gca,'fontsize',14);
axis tight;

%% Example 3 RF-Bode

clc; clear; close all;
z = 0.2;                    % damping ratio 
wrwn = sqrt(1-2*z^2);       % ratio of res freq. to natural freq. 
disp(wrwn); 

RF = tf([2*z 1],[1 2*z 1]); % take xfer function 
w = logspace(-1,2,1e4);     % driving freq range

% Get data for Bode 
[M, phi] = bode(RF,w); grid on;
M = M(:);   phi = phi(:);

% Plot Rel. Freq. Bode
subplot(2,1,1);
title('Bode Diagram','fontsize',14);
semilogx(w,20*log10(M),'b-');
grid on;
ylabel('TR [dB]','fontsize',14);
set(gca,'fontsize',14);
axis tight;

subplot(2,1,2);
semilogx(w,phi,'b-');
grid on;
ylabel('\phi [\circ]','fontsize',14);
xlabel('\omega/\omega_n','fontsize',14);
set(gca,'fontsize',14);
axis tight;

%% Example 3 Simulate
close all; clc; 

sim('l14_3.slx');
subplot(2,1,1);
plot(simout.time,simout.data(:,1),'k-');
grid on; axis tight;
title('Input','fontsize',14);
set(gca,'fontsize',14);
subplot(2,1,2);
plot(simout.time,simout.data(:,2),'k-');
grid on; axis tight;
title('Output','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

%% Example 4
clc; clear; close all;

RF1 = tf([1 0 0],[1 1 1]);
w = logspace(-2,2,2e4);     % driving freq range
bode(RF1,w); grid on;
xlabel('\omega/\omega_n');

RF2 = tf([1 0 0],[1 1.5 1]);
RF3 = tf([1 0 0],[1 0.75 1]);

bode([RF2 RF3],w); grid on;
xlabel('\omega/\omega_n');

%% Example 5: MDOF
clc; clear; close all;
k1 = 20; k2 = 30; 
b1 = 5; m1 = 10; m2 = 15;

A = [0 1 0 0
    -k1/m1 -b1/m1 k1/m1 b1/m1
    0 0 0 1
    k1/m2 b1/m2 -(k1+k2)/m2 -b1/m2];
lambda = eig(A);
disp(lambda);               % poles of system (2 CC pairs)
wn = abs(lambda([1,3]));    % natural freqs of system [rad/s]
disp(wn); 
% z = real(lambda([1,3]))./(-wn);

B = [0 1/m1 0 0]';
C = [1 0 0 0
    0 0 1 0];
D = [0 0]';

[b,a]=ss2tf(A,B,C,D);

pols = roots(a);    % find poles of TF
disp(abs(pols))     % natural freqs. of system 

TFz = tf(b(1,:),a);
TFy = tf(b(2,:),a);

clf;
bode(TFz,1e-1:0.01:1e1); grid on; hold on;
bode(TFy,1e-1:0.01:1e1);
legend('Z','Y');

Mz = 10^(-47.3/20); disp(Mz);
My = 10^(-67.5/20); disp(My);

[b,a]= ss2tf(A,B,C,D);
disp(b); disp(a);

Gz = tf(b(1,:),a);
Gy = tf(b(2,:),a);

%% Plot Example 5: MDOF
sim('l14_5.slx');   % run simulink model

figure(2); clf;     % full response 

tims = simout.time;
subplot(3,1,1);
plot(tims,simout.data(:,1),'k-');
grid on; axis tight;
title('u','fontsize',14);
set(gca,'fontsize',14);
subplot(3,1,2);
plot(tims,simout.data(:,2),'k-');
grid on; axis tight;
title('z','fontsize',14);
set(gca,'fontsize',14);
subplot(3,1,3);
plot(tims,simout.data(:,3),'k-');
grid on; axis tight;
title('y','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

figure(3); clf;     % zoomed to SS

sp = 3.8e6;
tims = simout.time(sp:end);
subplot(3,1,1);
plot(tims,simout.data(sp:end,1),'k-');
grid on; axis tight;
title('u','fontsize',14);
set(gca,'fontsize',14);
subplot(3,1,2);
plot(tims,simout.data(sp:end,2),'k-');
grid on; axis tight;
title('z','fontsize',14);
set(gca,'fontsize',14);
subplot(3,1,3);
plot(tims,simout.data(sp:end,3),'k-');
grid on; axis tight;
title('y','fontsize',14);
xlabel('Time [s]','fontsize',14);
set(gca,'fontsize',14);

%% Example 6 Node
clc; clear; close all;
k1 = 4e3; k2 = 3e3; 
b2 = 345; m = 20;

TR = tf([b2*(k1+k2) k1*k2],[b2*m k2*m b2*(k1+k2) k1*k2]);

clf;
bode(TR,1e0:0.1:1e2); grid on;
figure(2);
M = 10^(-11.2/20); disp(M);

den = [b2*m k2*m b2*(k1+k2) k1*k2];
p = roots(den); disp(p);
sqrt(abs(p))

tim = linspace(0,3,1e5);
u = 600*sin(40*tim);
lsim(TR,u,tim); grid on;
legend('F');

disp(10^(-11.2/20)); 

%% Example 7: LP Filter
clc; clear; close all;
wc = 5*2*pi; % cutoff frequency 

figure(1);
RF1 = tf(1,[1/(wc) 1]);
bode(RF1,1e-1:0.1:1e3); grid on;

figure(2);
t = linspace(0,15,1e4);
u = 5*sin(0.5*t)+0.6*sin(100*t);    % input signal 

y = lsim(RF1,u,t);                  % run input through filter block "RF1" for output signal "y"

plot(t,u,'b-',t,y,'r-');
xlabel('Time [s]','fontsize',14);
ylabel('Voltage [V]','fontsize',14);
legend('Unfiltered','Filtered','location','best');
set(gca,'fontsize',14);
grid on;

%% Example 8: HP Filter
clc; clear; close all;
wc = 5*2*pi; % cutoff frequency 

figure(1);
RF1 = tf([1/wc 0],[1/(wc) 1]);
bode(RF1,1e-1:0.1:1e3); grid on;

figure(2);
t = linspace(0,15,1e4);
u = 5*sin(0.5*t)+0.6*sin(100*t);

y = lsim(RF1,u,t);                  % run input through filter block "RF1" for output signal "y"

plot(t,u,'b-',t,y,'r-');
xlabel('Time [s]','fontsize',14);
ylabel('Voltage [V]','fontsize',14);
legend('Unfiltered','Filtered','location','best');
set(gca,'fontsize',14);
grid on;

%% Example 9: nth-order Passive LP Filter
clc; clear; close all;
wc = 5*2*pi; % cutoff frequency 

figure(1);
RF1 = tf(1,[1/(wc) 1]);
bode(RF1,1e-1:0.1:1e3); grid on;
hold on;
bode(RF1^2,1e-1:0.1:1e3); grid on;
bode(RF1^3,1e-1:0.1:1e3); grid on;
legend('1^{st}-Order','2^{nd}-Order','3^{rd}-Order','location','southwest');

figure(2);
t = linspace(0,15,1e4);
u = 5*sin(0.5*t)+0.6*sin(100*t);    % input 

y = lsim(RF1,u,t);                  % run input through filter block "RF1" for output signal "y"
y2 = lsim(RF1^2,u,t);               % run input through filter block "RF2" for output signal "y2"
y3 = lsim(RF1^3,u,t);               % run input through filter block "RF3" for output signal "y3"

plot(t,u,'b-',t,y,'r-',t,y2,'g-',t,y3,'c-');
xlabel('Time [s]','fontsize',14);
ylabel('Voltage [V]','fontsize',14);
legend('Unfiltered','1^{st}-Order','2^{nd}-Order','3^{rd}-Order','location','best');
set(gca,'fontsize',14);
grid on;
