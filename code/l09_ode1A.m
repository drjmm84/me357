%% State-space model for L09.1

function [xd] = l09_ode1A(t,x)

m = 4; % [kg]
k = 15; % [N/m]
b = 20; % [N-s/m]

% define input
% u = @(t) heaviside(t-3);  % delayed step function
u = @(t) 0;                 % no external input 

A = [ 
    0 1
    -k/m, -b/m
    ];

B = [0; 1/m];

xd = A*x + B*u(t);

return
