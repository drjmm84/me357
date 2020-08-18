%{
Runge-Kutta Example
%}

%% 1st-Order System Example  

clear; close all; clc;
clf; hold on;
% define derivative 
dx = @(x,t) -3*x;

st = 6:15:256;
for j = 1:length(st)
    t = linspace(0,2,st(j));
    dt = mean(diff(t));
    X = zeros(length(t),1);
    % IC
    X(1) = 2;
    
    for k=2:length(X)
        DX = dx(X(k-1),t(k-1));
        X(k) = X(k-1) + DX*dt;
    end
    
    plot(t,X);
    %     legend(num2str(dt));
end

plot(t,2*exp(-3*t),'k--');


grid on;
xlabel('t','fontsize',14);
ylabel('y(t)','fontsize',14);
set(gca,'fontsize',14);

% legend(num2str(st(1)),num2str(st(2)),num2str(st(3)),num2str(st(4)),...
%     num2str(st(5)),num2str(st(6)),'Analytic','location','northeast');

%% 2nd-Order System Example

clear; close all; clc;
clf; hold on;
% define derivatives 
dx1 = @(x,t) x(2);
dx2 = @(x,t) -2*x(1) - 1*x(2);

st = 30:60:800;
for j = 1:length(st)
    t = linspace(0,12,st(j));
    dt = mean(diff(t));
    X = zeros(length(t),2);
    % ICs
    X(1,1) = 2;
    X(1,2) = 0;
    for k=2:length(X)
        DX1 = dx1(X(k-1,:),t(k-1));
        X(k,1) = X(k-1,1) + DX1*dt;
        DX2 = dx2(X(k-1,:),t(k-1));
        X(k,2) = X(k-1,2) + DX2*dt;
    end
    
    plot(t,X(:,1));
    %     legend(num2str(dt));
end

% A = [0 1; -1 -2];
% l = eig(A);
plot(t,exp(-0.5*t).*(2.* cos(1.32288*t) + 0.755929*sin(1.32288.* t)),'k--')

grid on;
xlabel('t','fontsize',14);
ylabel('y(t)','fontsize',14);
set(gca,'fontsize',14);

% legend(num2str(st(1)),num2str(st(2)),num2str(st(3)),num2str(st(4)),...
%     num2str(st(5)),num2str(st(6)),num2str(st(7)),num2str(st(8)),...
%     num2str(st(9)),'location','northeast');