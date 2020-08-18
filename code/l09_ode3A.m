%% SS model for L09.3

function [xd] = l09_ode3A(t,x)

m1 = 5; m2 = 2; m3 = 3; % [kg]
k1=3; k2=2; k3=5;       % [N/m]
c1=6;                   % [N-s/m]
c2=c1; c3=c1; %#ok<NASGU>

% Define inputs 

f1 = @(t) heaviside(t-10);
% f1 = @(t) sin(t);
% f1 = @(t) 0;

% f2 = @(t) heaviside(t-5);
f2 = @(t) 0;
% f2 = @(t) sin(2*t);

% f3 = @(t) 5*heaviside(t);
f3 = @(t) 0;

u = [f1(t) f2(t) f3(t)]';

A = [ 
    0 1 0 0 0 0 0
    -k1/m1 -c1/m1 k1/m1 c1/m1 0 0 0
    0 0 0 1 0 0 0
    k1/m2 c1/m2 -(k1+k2)/m2 -c1/m2 k2/m2 0 0
    0 0 k2/c2 0 -k2/c2 0 1
    0 0 0 0 0 0 1
    0 0 k2/m3 0 -k2/m3 -k3/m3 0
    ];

B = [ 
    0 0 0
    1/m1 0 0
    0 0 0
    0 1/m2 0
    0 0 0
    0 0 0
    0 0 1/m3
    ];
    
xd = A*x + B*u;

return
