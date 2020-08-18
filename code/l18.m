%% Lecture 18

%% Find zeros for design 1 
clc; clear;

% use atan2 
fzero(@(b) atan2(-2.11+b,0.45)+atan2(-2.11-b,0.45)+7*pi/180,1)              % find imaginary component, b
% use angle 
fzero(@(b) angle(-1.15-2.11j + 1.6 + b*1j) + angle(-1.15-2.11j + 1.6 - b*1j) + 7*pi/180, 1.6)  % find imaginary component, b

%% Find zeros for design 2
clc; clear;

% fzero(@(b) atan2(-2.11+b,4.85)+atan2(-2.11-b,4.85)+7*pi/180,7)  % find imaginary component, b

fzero(@(b) angle(-1.15-2.11j + 6 + b*1j) + angle(-1.15-2.11j + 6 - b*1j) + 7*pi/180, 6)  % find imaginary component, b

%% PID
clc; clear; close all;

design =  2; % design 1 or 2
figure(1)
G0 = tf(1,[1 2]);
% choose design
if design == 1
    % root locus
    G1 = tf([1 3.2 22.3],[1 2 0]);
    rlocus(G1,0.01:0.001:1); grid on; hold on;
    axis equal;
    axis([-3,0,-5,5]);
    set(gca,'fontsize',12);
    % gains from root locus
    Kp = 1.09;
    Kd = 0.34;
    Ki = 7.57;
    GC = tf([Kd Kp Ki],[1e-6 1 0]);
    GCL = (GC*G0)/(1+GC*G0);
elseif design == 2
    % root locus
    G1 = tf([1 12 182],[1 2 0]);
    rlocus(G1,0.01:0.001:4); grid on; hold on;
    axis equal;
    axis([-3,0,-2.5,2.5]);
    set(gca,'fontsize',12);
    % gains from root locus
    Kp = 0.372;
    Kd = 0.031;
    Ki = 5.65;
    GC = tf([Kd Kp Ki],[1e-6 1 0]);
    GCL = (GC*G0)/(1+GC*G0);
end

% run simulation
sim('l18_ex1');
% get data from struct
tout = simout.time;
ref = simout.data(:,1);
eo = simout.data(:,3);
ei = simout.data(:,2);
err = simout.data(:,4);
% Plot results 
figure(2); clf;
subplot(2,1,1);
plot(tout,eo,'b-'); axis([0 max(tout) 0 6]);
ylabel('E_o [V]','fontsize',12);
xticks(0:1:max(tout));
yticks(-2:1:6);
set(gca,'fontsize',12);
grid on;
% plot input
subplot(2,1,2);
plot(tout,ei,'b-'); axis tight;
ylabel('E_i [V]','fontsize',12);
set(gca,'fontsize',12);
grid on;
% plot error
% subplot(3,1,2);
% plot(tout,err,'b-'); axis tight;
xlabel('Time [s]','fontsize',12);
% ylabel('Error','fontsize',12);
xticks(0:1:max(tout));
% yticks(-2:1:6);
% set(gca,'fontsize',12);
% grid on;
% find max overshot
[Mp] = max((eo-eo(end))/eo(end));
% find rise time (10%-90%)
indR = find(eo/eo(end)>0.9,1);  % time at 90 of SS%
indL = find(eo/eo(end)>0.1,1);  % time at 10 of SS%
tr = tout(indR)-tout(indL);     % rise time
disp([100*Mp tr])

return

%% TF Compare
% the CL TFs do not account for the saturation block in Simulink
close all; clc; clear; %#ok<*UNRCH>
TF1 = tf([0.253734 0.813433 5.64925],[1 2.30597 5.64925]);
TF2 = tf([0.00673434 0.322811 5.69141],[1 2.30934 5.69141]);
step(5*TF1);
hold on;
step(5*TF2);
legend('Design 1','Design 2','location','southeast');
