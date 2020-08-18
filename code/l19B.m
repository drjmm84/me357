%% L19 Nonlinear Dynamics

%% Lorenz System
% traces trajectories of 4 "close" ICs in Lorenz System
close all; clear; clc;
% System Parameters
rho = 28;
sigma = 10;
beta = 8/3;
T = 0:0.01:50;

% Simulation (Perturb ICs)
% 
% [~,X]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1 0 1]);
% % three *slightly* different ICs from the unperturbed system
% [~,X1]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1+1e-4*randn 0 1]);
% [~,X2]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1+1e-4*randn 0 1]);
% [~,X3]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1+1e-4*randn 0 1]);

% Simulation (Perturb rho)

[~,X]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1 0 1]);
% three *slightly* different ICs from the unperturbed system
[~,X1]= ode45(@(t,x) lorenz(t,x,rho*(1+1e-4*randn),sigma,beta),T,[1 0 1]);
[~,X2]= ode45(@(t,x) lorenz(t,x,rho*(1+1e-4*randn),sigma,beta),T,[1 0 1]);
[~,X3]= ode45(@(t,x) lorenz(t,x,rho*(1+1e-4*randn),sigma,beta),T,[1 0 1]);

% Initialize Movie
%vidObj = VideoWriter('ME357_16_A.mp4','MPEG-4');
%open(vidObj);

% Plotting
figure(1); clf;
for k=1:10:length(T)
    clf;
    plot3(X(1:k,1),X(1:k,2),X(1:k,3),'k-'); hold on;
    plot3(X(k,1),X(k,2),X(k,3),'go','MarkerSize',15,'MarkerFaceColor','g');
    plot3(X1(k,1),X1(k,2),X1(k,3),'md','MarkerSize',15,'MarkerFaceColor','m');
    plot3(X2(k,1),X2(k,2),X2(k,3),'cs','MarkerSize',15,'MarkerFaceColor','c');
    plot3(X3(k,1),X3(k,2),X3(k,3),'rp','MarkerSize',15,'MarkerFaceColor','r');
    % Adjust Visuals
    axis([floor(min(X(:,1))) ceil(max(X(:,1))) floor(min(X(:,2))) ...
        ceil(max(X(:,2))) floor(min(X(:,3))) ceil(max(X(:,3)))]);
    xlabel('x');    ylabel('y');    zlabel('z');
    grid on;
    set(gca,'fontsize',14);
    title(['t = ' num2str(T(k))]);
    drawnow;
    % Create Movie
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
    
end

figure(2); clf;
subplot(3,1,1);
plot(T,X(:,1),'k-',T,X1(:,1),'b-',T,X2(:,1),'r-',T,X3(:,1),'g-');
ylabel('x');
set(gca,'fontsize',14); grid on;
subplot(3,1,2);
plot(T,X(:,2),'k-',T,X1(:,2),'b-',T,X2(:,2),'r-',T,X3(:,2),'g-');
ylabel('y');
set(gca,'fontsize',14); grid on;
subplot(3,1,3);
plot(T,X(:,3),'k-',T,X1(:,3),'b-',T,X2(:,3),'r-',T,X3(:,3),'g-');
xlabel('Time'); ylabel('z');
set(gca,'fontsize',14); grid on;
%close(vidObj);

%% PDF of X value *at* T = 30
close all; clc;

rho = 28;
sigma = 10;
beta = 8/3;

x30 = zeros(1e5,1); % number of samples to use
x0 = x30;

% Monte Carlo Simulation of random start point *near* [1 0 1]
parfor k = 1:length(x30)
    [~,X]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),[0:30],[1+1e-4*randn 0 1]);
    x0(k) = X(1,1);     % save initial value 
    x30(k) = X(end,1);  % find x-value at 30 seconds
end
clear X* T;

figure(1); clf;
% subplot(2,1,1);
histogram(x0,'BinMethod','Scott','Normalization','Probability');
xlabel('$x(0)$','Interpreter','latex'); grid on;
set(gca,'fontsize',14);

disp(mean(x0));
disp(std(x0));

figure(2);
% subplot(2,1,2);
histogram(x30,'BinMethod','Scott','Normalization','Probability');
xlabel('$x(30)$','Interpreter','latex'); grid on;
set(gca,'fontsize',14);

disp(mean(x30));
disp(std(x30));

figure(3);
semilogx([1:length(x30)],cumsum(x30)./[1:length(x30)]'); %#ok<*NBRAK>
xlabel('$N$','Interpreter','latex');
ylabel('$\bar{x}(30)$','Interpreter','latex'); grid on;
set(gca,'fontsize',12);
axis tight;

%% PDF of X value *at* 30 seconds for SMD
close all; clc;

x30 = zeros(1e5,1);
x0 = x30;
% Monte Carlo Simulation of random start point *near* x0 = [1 0], m = 1, b = 1, k = 3
parfor k = 1:length(x30)
    [~,X]= ode45(@(t,y) l10_lin(t,y,1+1e-3*randn,1+1e-3*randn,3+1e-3*randn),0:60,[1+1e-3*randn 1e-3*randn]);
    x0(k) = X(1,1);
    x30(k) = X(end,1);
end
clear X* T;

figure(1); clf;
% subplot(2,1,1);
histogram(x0,'BinMethod','Scott','Normalization','Probability');
xlabel('$x(0)$','Interpreter','latex'); grid on;
set(gca,'fontsize',14);

disp(mean(x0));
disp(std(x0));

figure(2); clf;
% subplot(2,1,2);
histogram(x30,'BinMethod','Scott','Normalization','Probability');
xlabel('$x(30)$','Interpreter','latex'); grid on;
set(gca,'fontsize',14);

disp(mean(x30));
disp(std(x30));

figure(3);
semilogx([1:length(x30)],cumsum(x30)./[1:length(x30)]'); %#ok<*NBRAK>
xlabel('$N$','Interpreter','latex');
ylabel('$\bar{x}(30)$','Interpreter','latex'); grid on;
set(gca,'fontsize',12);
axis tight;

%% 95% Estimate X value *over* T = 45 
close all; clc;

xdata = zeros(2e3,2e3);
T = linspace(0,45,2e3);
parfor k=1:size(xdata,1)
    [~,X]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1+1e-4*randn 0 1]);
    xdata(k,:) = X(:,1);
end
clear X*;

figure(1); clf;
% xm = mean(xdata,1);
sl = quantile(xdata,0.025,1);
su = quantile(xdata,0.975,1);

plot(T,su,'r-',T,sl,'b-'); %#ok<*NBRAK>
xlabel('Time [s]');
ylabel('x'); grid on;
set(gca,'fontsize',14);
axis tight;
legend('97.5%','2.5%','location','best')

%% 95% Estimate X value *over* 20 seconds SMD
close all; clc;

xdata = zeros(5e3,2e3);
T = linspace(0,20,2e3);
parfor k=1:size(xdata,1)
    [~,X]= ode45(@(t,y) l10_lin(t,y,1+1e-3*randn,1+1e-3*randn,3+1e-3*randn),T,[1+1e-3*randn 1e-3*randn]);
    xdata(k,:) = X(:,1);
end
clear X*;

figure(1); clf;
% xm = mean(xdata,1);
sl = quantile(xdata,0.025,1);
su = quantile(xdata,0.975,1);

plot(T,su,'r-',T,sl,'b-'); %#ok<*NBRAK>
xlabel('Time [s]');
ylabel('x'); grid on;
set(gca,'fontsize',14);
axis tight;
legend('97.5%','2.5%','location','best')

%% 95% Estimate X value *over* 30 seconds re-estimate 
close all; clc;

rho = 28;
sigma = 10;
beta = 8/3;

T = linspace(0,30,2e3);
Tr = linspace(0,10,5e2);

xdata = zeros(1e3,length(T));
xdatar = zeros(1e3,length(Tr));

% Actual ICs (unknown at the start)
[~,xtst]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),0:20,[1.01 0 1]);
xtst = xtst(end,:);

parfor k=1:size(xdata,1)
    [~,X]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),T,[1+1e-4*randn 0 1]);
    xdata(k,:) = X(:,1);
    
    [~,Xr]= ode45(@(t,x) lorenz(t,x,rho,sigma,beta),Tr,xtst+[1e-4*randn 0 0]);
    xdatar(k,:) = Xr(:,1);
    
end
clear X*;

figure(1); clf;
% xm = mean(xdata,1);
sl = quantile(xdata,0.025,1);
su = quantile(xdata,0.975,1);

slr = quantile(xdatar,0.025,1);
sur = quantile(xdatar,0.975,1);

plot(T,su,'r-',T,sl,'b--'); %#ok<*NBRAK>

xlabel('Time');
ylabel('x'); grid on;
set(gca,'fontsize',14);
axis tight;
legend('97.5%','2.5%','location','best')

hold on; pause;
plot(20,xtst(1),'rs'); pause;
plot(Tr+20,sur,'y-',Tr+20,slr,'g--'); %#ok<*NBRAK>


