clear all; close all; clc

x1d = @(x1,x2) x2;
x2d = @(x1,x2) 4.*(1-x1.^2).*x2 - x1;

[X1, X2] = meshgrid(linspace(-6,6,25),linspace(-6,6,25));

X1d = x1d(X1,X2);
X2d = x2d(X1,X2);

quiver(X1,X2,X1d,X2d);
axis equal
axis tight
