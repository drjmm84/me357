%%

function [xd] = l10A_nonlin(~,x)

xd = [0 0]';

xd(1) = x(2);
xd(2) = -4*x(1)-0.5*abs( x(2) )*x(2);

return