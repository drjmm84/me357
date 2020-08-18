%% Lecture 9 Example 1B
clc; close all; clear;

m = 1; % kg
k = 16 ; % N/m
b = 10; % Ns/m

A = [ 0 1
    -k/m, -b/m];
B = [0; 1/m];
C = [1 0]; % pick out only y
D = 0;

[b,a] = ss2tf(A,B,C,D);

sys = tf(b,a)

disp(eig(A)); % eigenvalues of state matrix, A
