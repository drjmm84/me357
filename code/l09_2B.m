%% Lecture 9 Example 2B
clc; close all; clear;

m1 = 10; % kg
m2 = 1;
k1 = 100 ; % N/m
k2 = 30;
b1 = 20; % Ns/m

A = [ 0 1 0 0
    -k1/m1 -b1/m1 k1/m1 b1/m1
    0 0 0 1
    k1/m2 b1/m2 -(k1+k2)/m2 -b1/m2 ];

B = [0
    0
    0
    1/m2 ];

C = [1 0 0 0
    0 0 1 0]; % pick out z and y
D=[0
    0];

[b,a] = ss2tf(A,B,C,D);

sysZ = tf(b(1,:),a)
sysY = tf(b(2,:),a)

disp(eig(A));

