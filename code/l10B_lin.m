%%

function [xd] = l10B_lin(t,x)

global F; 

A = [0 1
    -1 0];

xd = A*x + [0;F]*cos(t);

return