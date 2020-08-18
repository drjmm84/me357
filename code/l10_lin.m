%% SS model for simple SMD

function [xd] = l10_lin(t,x,m,b,k)

A = [0 1
    -k/m -b/m];

xd = A*x;

return