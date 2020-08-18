%% Lecture 17
clear; close all; clc;

G0 = tf(1,[1 2 0]);
rlocus(G0); grid on; hold on;
scatter(-2,2,265,'rp');     % desired CL poles 
scatter(-2,-2,265,'rp');     % desired CL poles 
axis equal;
axis([-3,0,-2.5,2.5]);
set(gca,'fontsize',12);

%% Lead Lag
clc; clear; close all;

% choose design
design =  3;

if design == 1
    K=8;
    z1=2;
    p1=4;
        
elseif design == 2
    K = 16;
    z1= 3;
    p1= 8;
     
elseif design == 3
    K = 5.3;
    z1= 1;
    p1= 2.65;
        
end

% Root Locus
figure(1);
GRL1 = tf([1 z1],[1 2+p1 2*p1 0]); % general RL 
rlocus(GRL1,0.1:0.01:20); grid on; hold on;
axis equal;
axis([-3,0,-2.5,2.5]);
set(gca,'fontsize',12);

% Simulation 
sim('l17_ex1_LC');
tout = simout.time;
ref = simout.data(:,1);
eo = simout.data(:,3);
ei = simout.data(:,2);
err = simout.data(:,4);

% Plot
figure(2)
subplot(3,1,1); % response
plot(tout,eo,'b-'); axis tight;
ylabel('E_o','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,2); % input effort
plot(tout,ei,'b-'); axis tight;
ylabel('E_i','fontsize',12);
set(gca,'fontsize',12);
grid on;

subplot(3,1,3); % error 
plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
ylabel('Error','fontsize',12);
set(gca,'fontsize',12);
grid on;

% System Response 
[Mp, tp] = max((eo-eo(end))/eo(end)); % find maximum value
tp = tout(tp);
disp([100*Mp tp max(ei)])

%% Compare
clf; clc; close all;

GPD = tf([1 4],[1 2 0]);            % PD Controller RL TF
GRL1 = tf([1 2],[1 6 8 0]);         % LC Controller RL TFs
GRL2 = tf([1 3],[1 10 16 0]);
GRL3 = tf([1 1],[1 4.65 5.3 0]);

GPD = feedback(2*GPD,tf(1,1));      % CL TF for PD
G1 = feedback(8*GRL1,tf(1,1));      % CL TF for LC
G2 = feedback(16*GRL2,tf(1,1));
G3 = feedback(5.3*GRL3,tf(1,1));

% Step Response 
figure(1); clf;
step(5*G1,'k-'); hold on; grid on;
step(5*G2,'r-');
step(5*G3,'g-');
step(5*GPD,'b-');
title(' ');
legend('G_1','G_2','G_3','G_{PD}','location','southeast');
set(gca,'fontsize',14);

% Bode Plots 
figure(2); clf;
bode(G1,'k-'); hold on; grid on;
bode(G2,'r-');
bode(G3,'g-');
bode(GPD,'b-');
legend('G_1','G_2','G_3','G_{PD}','location','southwest');
set(gca,'fontsize',14);
title(' ');

%% Zero and Pole Choices 
clf; clc; close all;
zc = linspace(0,3.0,25);
pcl = zeros(1,length(zc));
for k=1:length(zc)
    pcl(k) = fzero(@(pc) -atan2(-2,-2+pc) + atan2(-2,-2+zc(k)) + 45*pi/180,8);
end
plot(zc,pcl,'k*-');
xlabel('z_c','fontsize',14); ylabel('p_c','fontsize',14);
set(gca,'fontsize',14);
axis tight;
grid on;

%% Compare CL
clc; clear; close all;
t = linspace(0,10,2e3);

TF1 = tf(8,[1 4 8]);
TF2 = tf([16 48],[1 10 32 48]);
TF3 = tf([5.3 5.3],[1 4.65 10.6 5.3]);

step(5*TF1,t);
hold on;
step(5*TF2,t);
step(5*TF3,t);

legend('Design 1','Design 2','Design 3','location','southeast');
set(gca,'fontsize',14);
grid on;
