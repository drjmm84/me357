%%

function [xd] = l10_nonlin(~,x)

xd = [0 0]';

xd(1) = x(2);
xd(2) = -3*x(1)-5*x(1)^3 - 1*x(2);

return