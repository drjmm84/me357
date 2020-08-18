%% State-space model for VDP 

function [xd] = l16_vdp(~,x,mu)

xd=[0;0];

xd(1) = x(2);
xd(2) = mu*(1-x(1)^2)*x(2)-x(1);

return