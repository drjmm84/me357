%%

function [xd] = l10B_nonlin(t,x)

global F;

xd = [0 0]';

xd(1) = x(2);
xd(2) = 0.5*x(1)-0.5*x(1)^3 + F*cos(t);

return