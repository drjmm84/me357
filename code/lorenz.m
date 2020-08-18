function [xd] = lorenz(~,x,rho, sigma, beta)

xd = [0;0;0];                   % preallocate size and shape of \dot{x_} 

xd(1) = sigma*(x(2)-x(1));      % \dot{x}
xd(2) = x(1)*(rho-x(3))-x(2);   % \dot{y}
xd(3) = x(1)*x(2)-beta*x(3);    % \dot{z}

return;