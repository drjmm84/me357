%%

function [xd] = l10A_lin(~,x)

A = [0 1
    -4 0];

xd = A*x;

return