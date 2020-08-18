%% L19 Nonlinear Dynamics

%% Linear Example (SMD)
clear; close all; clc;

% ICs
x0 = 1;
xd0 = 0;
% Parameters
k = 3;
b = 1;
m = 1;

Ti = 0:0.1:12;  % simulate for 12 seconds

[T,X] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0 xd0]);

% Plot Time Trace
figure(1); clf;
for k=1:length(T)
    
    subplot(2,1,1);
    plot(T(1:k),X(1:k,1),'k-');
    axis([0 max(T) floor(min(X(:,1))) ceil(max(X(:,1))) ]);
    ylabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Position','fontsize',12);
    
    subplot(2,1,2);
    plot(T(1:k),X(1:k,2),'k-');
    axis([0 max(T) floor(min(X(:,2))) ceil(max(X(:,2))) ]);
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('time [s]','fontsize',12);
    set(gca,'fontsize',12);
    grid on;
    title('Velocity','fontsize',12);
    drawnow
    
end

pause;

% Plot Phase Plane
figure(2); clf;
for k=1:length(T)
    
    plot(X(1:k,1),X(1:k,2),'k-');
    axis square;
    axis([floor(min(X(:,1))) ceil(max(X(:,1))) ...
        floor(min(X(:,2))) ceil(max(X(:,2)))]);
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Phase Plane','fontsize',12);
    drawnow
    
end
hold on;
plot(X(end,1),X(end,2),'rx');

%% Linear System with IC perturbations

% Plot Phase Plane

% ICs
x0 = 1;
xd0 = 0;
% Parameters
k = 3;
b = 1;
m = 1;

Ti = 0:0.1:12;

% run perturbations off baseline IC
[T,X] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0 xd0]);   % baseline
[~,X1] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X2] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X3] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X4] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X5] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X6] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);
[~,X7] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0+0.075*randn xd0 + 0.075*randn]);

figure(3); clf;
for k=1:length(T)
    
    subplot(2,1,1);
    plot(T(1:k),X(1:k,1),'k-'); hold on;
    plot(T(1:k),X1(1:k,1),'r-');
    plot(T(1:k),X2(1:k,1),'g-');
    plot(T(1:k),X3(1:k,1),'c-');
    plot(T(1:k),X4(1:k,1),'m-');
    plot(T(1:k),X5(1:k,1),'b-');
    plot(T(1:k),X6(1:k,1),'y-');
    plot(T(1:k),X7(1:k,1),'k:');
    
    axis([0 max(T) floor(min(X(:,1))) 1.2*ceil(max(X(:,1))) ]);
    ylabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Position','fontsize',12);
    
    subplot(2,1,2);
    plot(T(1:k),X(1:k,2),'k-'); hold on;
    plot(T(1:k),X1(1:k,2),'r-');
    plot(T(1:k),X2(1:k,2),'g-');
    plot(T(1:k),X3(1:k,2),'c-');
    plot(T(1:k),X4(1:k,2),'m-');
    plot(T(1:k),X5(1:k,2),'b-');
    plot(T(1:k),X6(1:k,2),'y-');
    plot(T(1:k),X7(1:k,2),'k:');
    
    axis([0 max(T) floor(min(X(:,2))) 1.2*ceil(max(X(:,2))) ]);
    
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('time [s]','fontsize',12);
    set(gca,'fontsize',12);
    grid on;
    title('Velocity','fontsize',12);
    drawnow
    
end

pause;

figure(4); clf;
for k=1:length(T)
    
    plot(X(1:k,1),X(1:k,2),'k-'); hold on;
    plot(X1(1:k,1),X1(1:k,2),'r-');
    plot(X2(1:k,1),X2(1:k,2),'g-');
    plot(X3(1:k,1),X3(1:k,2),'c-');
    plot(X4(1:k,1),X4(1:k,2),'m-');
    plot(X5(1:k,1),X5(1:k,2),'b-');
    plot(X6(1:k,1),X6(1:k,2),'y-');
    plot(X7(1:k,1),X7(1:k,2),'k:');
    
    axis square;
    axis([floor(min(X(:,1))) 1.2*ceil(max(X(:,1))) ...
        floor(min(X(:,2))) 1.2*ceil(max(X(:,2)))]);
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Phase Plane','fontsize',12);
    drawnow
    
end

%% Linear System with parameter perturbations

% Plot Phase Plane

% ICs
x0 = 1;
xd0 = 0;
% Parameters
k = 3;
b = 1;
m = 1;

Ti = 0:0.1:12;

% run perturbations off baseline parameters
[T,X] = ode45(@(t,y) l10_lin(t,y,m,b,k),Ti,[x0 xd0]);   % baseline
[~,X1] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X2] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X3] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X4] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X5] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X6] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);
[~,X7] = ode45(@(t,y) l10_lin(t,y,m+0.2*randn,b+0.2*randn,k+0.2*randn),Ti,[x0 xd0]);

figure(5); clf;
for k=1:length(T)
    
    subplot(2,1,1);
    plot(T(1:k),X(1:k,1),'k-'); hold on;
    plot(T(1:k),X1(1:k,1),'r-');
    plot(T(1:k),X2(1:k,1),'g-');
    plot(T(1:k),X3(1:k,1),'c-');
    plot(T(1:k),X4(1:k,1),'m-');
    plot(T(1:k),X5(1:k,1),'b-');
    plot(T(1:k),X6(1:k,1),'y-');
    plot(T(1:k),X7(1:k,1),'k:');
    
    axis([0 max(T) floor(min(X(:,1))) 1.2*ceil(max(X(:,1))) ]);
    ylabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Position','fontsize',12);
    
    subplot(2,1,2);
    plot(T(1:k),X(1:k,2),'k-'); hold on;
    plot(T(1:k),X1(1:k,2),'r-');
    plot(T(1:k),X2(1:k,2),'g-');
    plot(T(1:k),X3(1:k,2),'c-');
    plot(T(1:k),X4(1:k,2),'m-');
    plot(T(1:k),X5(1:k,2),'b-');
    plot(T(1:k),X6(1:k,2),'y-');
    plot(T(1:k),X7(1:k,2),'k:');
    
    axis([0 max(T) floor(min(X(:,2))) 1.2*ceil(max(X(:,2))) ]);
    
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('time [s]','fontsize',12);
    set(gca,'fontsize',12);
    grid on;
    title('Velocity','fontsize',12);
    drawnow
    
end

pause;

figure(6); clf;
for k=1:length(T)
   
    plot(X(1:k,1),X(1:k,2),'k-'); hold on;
    plot(X1(1:k,1),X1(1:k,2),'r-');
    plot(X2(1:k,1),X2(1:k,2),'g-');
    plot(X3(1:k,1),X3(1:k,2),'c-');
    plot(X4(1:k,1),X4(1:k,2),'m-');
    plot(X5(1:k,1),X5(1:k,2),'b-');
    plot(X6(1:k,1),X6(1:k,2),'y-');
    plot(X7(1:k,1),X7(1:k,2),'k:');
    
    axis square;
    axis([floor(min(X(:,1))) 1.2*ceil(max(X(:,1))) ...
        floor(min(X(:,2))) 1.2*ceil(max(X(:,2)))]);
    ylabel('$\dot{x}(t)$','fontsize',14,'interpreter','latex');
    xlabel('${x}(t)$','fontsize',14,'interpreter','latex');
    set(gca,'fontsize',12);
    grid on;
    title('Phase Plane','fontsize',12);
    drawnow
    
end
