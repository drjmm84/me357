%% Lecture 15 driver file

%% OL Step
clear; close all; clc;      %#ok<*CLALL>

G = tf(1,[2.5 1]);
t = linspace(0,20,1e3);     % run for 20 sec

% step input
u1 = 5*ones(length(t),1);
y1 = lsim(G,u1,t);
figure(1);
plot(t,u1,'b',t,y1,'r');
legend('E_i','E_o','location','southeast');
set(gca,'fontsize',12);
ylabel('Voltage','fontsize',12);
xlabel('Time [s]','fontsize',12);
grid on;

% ramp input
u2 = min(t,5);
y2 = lsim(G,u2,t);
figure(2);
plot(t,u2,'b',t,y2,'r');
legend('E_i','E_o','location','southeast');
set(gca,'fontsize',12);
ylabel('Voltage','fontsize',12);
xlabel('Time [s]','fontsize',12);
grid on;

%% P Controller
clear; close all; clc; 

Kp = 10;                 % Gain on proportion 
sim('l15_ex1_P');
tout = simout.time;
ref = simout.data(:,1); % Reference signal
eo = simout.data(:,2);  % Output Voltage
ei = simout.data(:,3);  % Input Voltage 
err = simout.data(:,4); % Error 

% Plot
subplot(4,1,1);
plot(tout,ref,'b-'); axis tight;
ylabel('Ref','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,2);
plot(tout,eo,'b-'); axis tight;
ylabel('E_o','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,3);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,4);
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
grid on;

%% P Controller Comparison
clear; close all; clc;

Kp_list = [1 5 10];                 % Gain on proportion
for k=1:3
    Kp = Kp_list(k);
    sim('l15_ex1_P');
    data(k).tout = simout.time;     %#ok<*SAGROW>
    data(k).ref = simout.data(:,1);
    data(k).eo = simout.data(:,2);
    data(k).ei = simout.data(:,3);
    data(k).err = simout.data(:,4);
end

subplot(2,1,1);
plot(data(1).tout,data(1).ei,'k-',data(2).tout,data(2).ei,'r-',data(3).tout,data(3).ei,'g-');
grid on;
ylabel('E_i [V]','fontsize',12);
set(gca,'fontsize',12);
legend('K_P=1','K_P=5','K_P=10','location','northeast')

subplot(2,1,2);
plot(data(1).tout,data(1).eo,'k-',data(2).tout,data(2).eo,'r-',data(3).tout,data(3).eo,'g-');
grid on;
xlabel('Time [s]','fontsize',12);
ylabel('E_o [V]','fontsize',12);
set(gca,'fontsize',12);

%% I Controller
clear; close all; clc; 

Ki = 0.1; % Gain on Integral 
wn = sqrt(2*Ki/5); disp(wn);
zeta = sqrt(1/(10*Ki)); disp(zeta);
wd = wn*sqrt(1-zeta^2);
tp = pi/wd; disp(tp);

sim('l15_ex2_I');
tout = simout.time;
ref = simout.data(:,1);
eo = simout.data(:,2);
ei = simout.data(:,3);
err = simout.data(:,4);
% Plot
subplot(4,1,1);
plot(tout,ref,'b-'); axis tight;
ylabel('Ref','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,2);
plot(tout,eo,'b-'); axis tight;
ylabel('E_o','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,3);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(4,1,4);
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
grid on;

%% I Controller Comparison
clear; close all; clc;
Ki_list = [1 5 25]; % Gain on proportion
for k=1:3
    Ki = Ki_list(k);
    sim('l15_ex2_I');
    data(k).tout = simout.time;
    data(k).ref = simout.data(:,1);
    data(k).eo = simout.data(:,2);
    data(k).ei = simout.data(:,3);
    data(k).err = simout.data(:,4);
end

subplot(2,1,1);
plot(data(1).tout,data(1).ei,'k-',data(2).tout,data(2).ei,'r-',data(3).tout,data(3).ei,'g-');
grid on;
ylabel('E_i [V]','fontsize',12);
%legend('K_I=1','K_I=5','K_I=25','location','southeast')
set(gca,'fontsize',12);

subplot(2,1,2);
plot(data(1).tout,data(1).eo,'k-',data(2).tout,data(2).eo,'r-',data(3).tout,data(3).eo,'g-');
grid on;
xlabel('Time [s]','fontsize',12);
ylabel('E_o [V]','fontsize',12);
legend('K_I=1','K_I=5','K_I=25','location','southeast')
set(gca,'fontsize',12);
