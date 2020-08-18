%% Lecture 16 driver 

%% OL
clear; close all; clc;  %#ok<*CLALL>

figure(1);
roots([1 2 0])
G0 = tf(1,[1 2 0]);
step(G0)

% add P control 
figure(2);
K = 1;
roots([1 2 K])
GCL = tf(K,[1 2 K]);
step(GCL)

%% P Controller
clear; close all; clc; 

Kp = 1; %#ok<NASGU>
sim('l16_ex1_P');
tout = simout.time;
%ref = simout.data(:,1);
eo = simout.data(:,2);
ei = simout.data(:,3);
err = simout.data(:,4);
% Plot

subplot(3,1,1);
plot(tout,eo,'b-'); axis tight;
ylabel('E_o','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,2);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,3);
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
grid on;

%% P Controller Comparison
clear; close all; clc;
Kp_list = [1 5 10]; % Gain on proportion
for k=1:3
    Kp = Kp_list(k); %#ok<NASGU>
    sim('l16_ex1_P');
    data(k).tout = simout.time; %#ok<*SAGROW>
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

subplot(2,1,2);
plot(data(1).tout,data(1).eo,'k-',data(2).tout,data(2).eo,'r-',data(3).tout,data(3).eo,'g-');
grid on;
xlabel('Time [s]','fontsize',12);
ylabel('E_o [V]','fontsize',12);
legend('K_P=1','K_P=5','K_P=10','location','best')
set(gca,'fontsize',12);

%% P Controller Root Locus
close all; clear; clc;

G0 = tf(1,[1 2 0]);     % This is G' = G0 = GOL 

rlocus(G0,1:1:15);      % Use G' in RL plot
grid on;
set(gca,'fontsize',12);
axis([-2.5 0 -4 4])

%% P Controller
clear; close all; clc; 

Kp = 4; %#ok<NASGU>
sim('l16_ex1_P');
tout = simout.time;
%ref = simout.data(:,1);
eo = simout.data(:,2);
ei = simout.data(:,3);
err = simout.data(:,4);
% Plot

subplot(3,1,1);
plot(tout,eo,'b-'); axis tight;
ylabel('E_o','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,2);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,3);
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
grid on;

%% finding zero
clc;

% disp((-atan2(-1.73,-1+2) - atan2(-1.73,-1+0))*180/pi);
disp( (-angle(-1-1.73j+2) - angle(-1-1.73j+0))*180/pi );

% angle from OL poles
% disp((-atan2(-2,-2+2) - atan2(-2,-2+0))*180/pi);
disp( (-angle(-2-2j+2) -angle(-2-2j+0))*180/pi );

% check zc choice 
disp( (-angle(-2-2j+2) -angle(-2-2j+0) +angle(-2-2j+4))*180/pi );

%% PD Example
clc; close all; clear;
figure(1)
G0 = tf(1,[1 2 0]);

% Root Locus 
GRL = tf([1 4],[1 2 0]); % G'(s) = (s+z)*G0(s)
rlocus(GRL,0:0.01:12);
grid on;
set(gca,'fontsize',12);

Kd = 2; % derivative gain 
Kp = 8; % proportional gain

figure(2)
sim('l16_ex2_PD');
tout = simout.time;
ref = simout.data(:,1);
eo = simout.data(:,2);
ei = simout.data(:,3);
err = simout.data(:,4);
% Plot

subplot(3,1,1);
plot(tout,eo,'b-'); axis tight;
ylabel('E_o [V]','fontsize',12);
set(gca,'fontsize',12);
xticks(0:1:6);
yticks(-2:1:6);
grid on;

subplot(3,1,3);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,2);
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
xticks(0:1:6);
yticks(-2:1:6);
grid on;

[Mp, tr] = max((eo-eo(end))/eo(end));   % max overshoot and index 
tr = tout(tr);                          % peak time
disp([100*Mp tr])

%% Compare
clear; close all; clc;

TF1 = tf(8,[1 4 8]);
TF2 = tf([2 8],[1 4 8]);

step(5*TF1);
hold on;
step(5*TF2);
grid on;

legend('Baseline','Extra Zero','location','southeast');
