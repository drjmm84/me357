%% State-space model for SMD 

function [xd] = l16_smd(t,x,m,k,b)

u = @(t) 0;

A = [ 0 1
    -k/m, -b/m];

B = [0; 1/m];

xd = A*x + B *u(t);

return