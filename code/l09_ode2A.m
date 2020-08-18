%% State-space model for L09.2

function [xd] = l09_ode2A(t,x)

m1 = 5; m2 = 2; % [kg]
k1=3; k2=2;     % [N/m]
b1=6;           % [N-s/m]

% Define input 
% u = @(t) heaviside(t-10); % delayed step function
% u = @(t) sin(t);          % sin input 
u = @(t) 0;                 % no external input

A = [ 
    0 1 0 0
    -k1/m1 -b1/m1 k1/m1 b1/m1
    0 0 0 1
    k1/m2 b1/m2 -(k1+k2)/m2 -b1/m2 ];

B = [0
    0
    0
    1/m2 ];

xd = A*x + B*u(t);

return
