%% Lecture 7 Catheter Example
clc; clear; close all;

l = 3;          %#ok<*NASGU> % tube length [m]
D = 5e-3;       % tube diameter [m]
C = 2.25e-12;   % sensor capacitance in [m^3/Pa], 1 mLHg = 7.50061683e-9 m^3/Pa

% second parameters
% l = 1;    % [m]
% D = 8e-3; % [m]

mu = 3.5e-3;    % viscosity of blood [N-s/m^2] | https://en.wikipedia.org/wiki/Hemorheology#Normal_level
rho = 1125;     % density of blood [kg/m^3] | https://hypertextbook.com/facts/2004/MichaelShmukler.shtml

R = 128*mu*l/(pi*D^4);  % tube wall resistance
L = rho*l/(pi*(D/2)^2); % blood inertia

sim('l07_fluidlevel_LV_catheter');    % run simulink model with current l, D, C
%  produces "tout" and "simout"

data = simout.Data; % get data from simulink model

PLV = data(:,1);    % LVP actual signal
PS = data(:,2);     % signal from pressure sensor;

% Plot response
plot(tout,PLV,'k',tout,PS,'r');
grid on;
legend('P_{LV}','P_{sensor}','location','northeast');
xlabel('Time [s]','fontsize',14);
ylabel('Pressure [Pa]','fontsize',14);
set(gca,'fontsize',14);

% RMS error | https://en.wikipedia.org/wiki/Root-mean-square_deviation
err = norm(PLV-PS)/sqrt(length(PS));
disp(err);

opts = optimoptions('fmincon','OptimalityTolerance',1e-12,'ConstraintTolerance',1e-12);
p = fmincon(@(X) designcheck(X),[9 8e-3],[],[],[],[],[0.5 1e-4],[10 1e-2],[],opts)

function [err] = designcheck(X)

l = X(1);
D = X(2);

mu = 3.5e-3;    % viscosity of blood [N-s/m^2] | https://en.wikipedia.org/wiki/Hemorheology#Normal_level
rho = 1125;     % density of blood [kg/m^3] | https://hypertextbook.com/facts/2004/MichaelShmukler.shtml

R = 128*mu*l/(pi*D^4);  % tube wall resistance
L = rho*l/(pi*(D/2)^2); % blood inertia

sim('l07_fluidlevel_LV_catheter');    % run simulink model with current l, D, C
data = simout.Data; % get data from simulink model
PLV = data(:,1);    % LVP actual signal
PS = data(:,2);     % signal from pressure sensor;
err = norm(PLV-PS)/sqrt(length(PS));

end
