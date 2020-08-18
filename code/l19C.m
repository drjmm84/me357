%% L19 Nonlinear Dynamics

%% Logistic Map Examples
clc; clear; close all;

r= 3.2;

x = zeros(5e3,length(r));

x(1) = rand; % random start point on [0,1]
for JJ = 2:length(x)
    x(JJ) = r*x(JJ-1)*(1-x(JJ-1)); %#ok<*PFBNS>
end
% rang=1:40;
rang = length(x)-4e3+500:length(x); % trim to steady state values

plot(rang,x(rang),'ko','MarkerSize',6);
axis([min(rang) max(rang) 0 1])
xlabel('N'); grid on;
ylabel('x_N')
set(gca,'fontsize',14);

%% random plot
figure(1); clf;
xr = rand(3e3,1);
plot(xr,'ko','MarkerSize',6);
xlabel('N'); grid on;
ylabel('x_N')
set(gca,'fontsize',14);
figure(2); clf;
histogram(xr,'BinMethod','Sturges','Normalization','Probability');
set(gca,'fontsize',14); grid on;
figure(3); clf;
scatter(xr(1:end-1),xr(2:end),'k.');
xlabel('x_N'); grid on;
ylabel('x_{N+1}')
set(gca,'fontsize',14); grid on;

%% Histogram
figure(2); clf;
histogram(x(rang),'BinMethod','Sturges','Normalization','Probability');
set(gca,'fontsize',14); grid on;
figure(3); clf;
scatter(x(4999:end-1),x(5000:end),'k.');
xlabel('x_N'); grid on;
ylabel('x_{N+1}')
set(gca,'fontsize',14); grid on;

%% Logistic Map: Bifurcation Diagram

clc; clear; close all;
r = linspace(2.8,4,1444);
xf = zeros(4e2,length(r));

for KK = 1:length(r)        % cycle different r-values
    for GG=1:size(xf,1)  % repeat multiple ICs for each r
        x = zeros(1e3,1);
        x(1) = rand;        % initial value is random on [0,1]
        for JJ = 2:length(x)
            x(JJ) = r(KK)*x(JJ-1)*(1-x(JJ-1)); %#ok<*PFBNS>
        end
        xf(GG,KK) = x(end); % save last x-value only
    end
end
% plot diagram
plot(r,xf,'k.','MarkerSize',1e-6);
axis equal;
axis([min(r) max(r) 0 1]);
xlabel('r');
ylabel('x_{ss}')
set(gca,'fontsize',12); grid on;
%%
hist(x(rang),100,'k')

scatter(x(8e4:9e4), x(8e4+1:9e4+1),'ko')
xlabel('x_n'); ylabel('x_{n+1}');
grid on;

%% SMD Phase Space

clf; clear; clc;

m = 2; % kg
k = 12; % N/m
b = 25; % Ns/m
ft = 45;

A = [ 0 1
    -k/m, -b/m];
[~,L] = eig(A);
disp(L);

figure(1); hold on;
for k =1:30
    [T,X]= ode45(@(t,x) l16_smd(t,x,m,k,b),0:0.01:ft,-5+10*rand(2,1));
    scatter(X(:,1),X(:,2),3,T);
end
grid on;
set(gca,'fontsize',14); axis equal;
xlabel('$x$','interpreter','latex');
ylabel('$\dot{x}$','interpreter','latex');
colorbar; colormap('jet');

%% SMD
clf; clear; clc;

m = 4; % kg
k = 5; % N/m
b = 18; % Ns/m
ft = 5;

A = [ 0 1
    -k/m, -b/m];
[V,L] = eig(A) %#ok<*NOPTS>

[T,X]= ode45(@(t,x) l16_smd(t,x,m,k,b),0:0.01:ft,2*[randn randn]);
[~,X1]= ode45(@(t,x) l16_smd(t,x,m,k,b),0:0.01:ft,2*[randn randn]);
[~,X2]= ode45(@(t,x) l16_smd(t,x,m,k,b),0:0.01:ft,2*[randn randn]);
[~,X3]= ode45(@(t,x) l16_smd(t,x,m,k,b),0:0.01:ft,2*[randn randn]);

x=linspace(-5,5,50);
v =x;
[x,v] = meshgrid(x,v);
dx = v;
dv = -k/m*x - b/m*v;

for k=1:5:length(T)
    clf;
    
    plot(X(1:k,1),X(1:k,2),'g-'); hold on;
    plot(X(k,1),X(k,2),'go','MarkerSize',15,'MarkerFaceColor','g');
    
    plot(X1(1:k,1),X1(1:k,2),'m-'); hold on;
    plot(X1(k,1),X1(k,2),'md','MarkerSize',15,'MarkerFaceColor','m');
    plot(X2(1:k,1),X2(1:k,2),'c-'); hold on;
    plot(X2(k,1),X2(k,2),'cs','MarkerSize',15,'MarkerFaceColor','c');
    plot(X3(1:k,1),X3(1:k,2),'r-'); hold on;
    plot(X3(k,1),X3(k,2),'rp','MarkerSize',15,'MarkerFaceColor','r');
    
    quiver(x,v,dx,dv,'k');
    set(gca,'fontsize',14); axis equal;
    axis([-5 5 -5 5]);
    xlabel('x');    ylabel('v');
    grid on;
    title(['t = ' num2str(T(k))]);
    pause(1e-6);
end

% Stability Manifolds
% plot(0.1*V(:,1),100*V(:,1),'r--')
% plot(0.1*V(:,2),100*V(:,2),'b--')

%% Duffing Osc.
close all;
figure(1);
clf; clear; clc;

alpha = 0; %
beta = 1;
delta = 1;
ft = 15;

A0 = [ 0 1; alpha -delta]; % Jacobian at 0
A1 = [0 1; -2*alpha -delta]; % Jacobian at +/- sqrt(alpha/beta)

l0 = eig(A0);
l1 = eig(A1);

disp(sqrt(alpha/beta));
disp(l0);
disp(l1);

T = 0:0.01:ft;

[T,X]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X1]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X2]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X3]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X4]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X5]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);
[~,X6]= ode45(@(t,x) l16_duff(t,x,alpha,beta,delta),T,4*[randn randn]);

for k=1:10:length(T)
    clf;
    
    plot(X(1:k,1),X(1:k,2),'g-'); hold on;
    plot(X(k,1),X(k,2),'go','MarkerSize',15,'MarkerFaceColor','g');
    plot(X1(1:k,1),X1(1:k,2),'m-'); hold on;
    plot(X1(k,1),X1(k,2),'md','MarkerSize',15,'MarkerFaceColor','m');
    plot(X2(1:k,1),X2(1:k,2),'c-'); hold on;
    plot(X2(k,1),X2(k,2),'cs','MarkerSize',15,'MarkerFaceColor','c');
    plot(X3(1:k,1),X3(1:k,2),'r-'); hold on;
    plot(X3(k,1),X3(k,2),'rp','MarkerSize',15,'MarkerFaceColor','r');
    plot(X4(1:k,1),X4(1:k,2),'b-'); hold on;
    plot(X4(k,1),X4(k,2),'b*','MarkerSize',15,'MarkerFaceColor','b');
    plot(X5(1:k,1),X5(1:k,2),'k-'); hold on;
    plot(X5(k,1),X5(k,2),'kx','MarkerSize',15,'MarkerFaceColor','k');
    plot(X6(1:k,1),X6(1:k,2),'y-'); hold on;
    plot(X6(k,1),X6(k,2),'y^','MarkerSize',15,'MarkerFaceColor','y');
    %     quiver(x,v,dx,dv,'k');
    set(gca,'fontsize',14); axis tight;
    %     axis(10*[-1 1 -1 1]);
    xlabel('x_1');    ylabel('x_2');
    grid on;
    title(['t = ' num2str(T(k))]);
    pause(1e-6);
end

%% VDP Osc.
figure(1);
clf; clear; clc;

mu = 0.725; %
ft = 60;

[T,X]= ode45(@(t,x) l16_vdp(t,x,mu),0:0.01:ft,0.2*[randn randn]);
[~,X1]= ode45(@(t,x) l16_vdp(t,x,mu),0:0.01:ft,1*[randn randn]);
[~,X2]= ode45(@(t,x) l16_vdp(t,x,mu),0:0.01:ft,2*[randn randn]);
[~,X3]= ode45(@(t,x) l16_vdp(t,x,mu),0:0.01:ft,4*[randn randn]);

x=linspace(-10,10,50);
v =x;
[x,v] = meshgrid(x,v);
dx = v;
dv = -2*mu*x.*v.*x + v.*(1-x.^2)*mu;

for k=1:20:length(T)
    clf;
    
    plot(X(1:k,1),X(1:k,2),'g-'); hold on;
    plot(X(k,1),X(k,2),'go','MarkerSize',15,'MarkerFaceColor','g');
    
    plot(X1(1:k,1),X1(1:k,2),'m-'); hold on;
    plot(X1(k,1),X1(k,2),'md','MarkerSize',15,'MarkerFaceColor','m');
    plot(X2(1:k,1),X2(1:k,2),'c-'); hold on;
    plot(X2(k,1),X2(k,2),'cs','MarkerSize',15,'MarkerFaceColor','c');
    plot(X3(1:k,1),X3(1:k,2),'r-'); hold on;
    plot(X3(k,1),X3(k,2),'rp','MarkerSize',15,'MarkerFaceColor','r');
    
    %     quiver(x,v,dx,dv,'k');
    axis equal;
    axis(5*[-1 1 -1 1]);
    xlabel('x_1');    ylabel('x_2');
    grid on;
    title(['t = ' num2str(T(k))]);
    pause(1e-6);
end

%% VDP mu
clc; clear; close all;

mu = linspace(-4,4,1e3);

l = 0.5*(mu + sqrt(mu.^2-4));

plot3(mu, real(l),imag(l),'k-');
grid;
xlabel('\mu');
ylabel('Re(\lambda)');
zlabel('Im(\lambda)');

