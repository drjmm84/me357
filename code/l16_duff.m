%% State-space model for VDP 

function [xd] = l16_duff(~,x,alpha,beta,delta)

xd=[0;0];

xd(1) = x(2);
xd(2) = -delta*x(2)+alpha*x(1)-beta*x(1)^3;

return