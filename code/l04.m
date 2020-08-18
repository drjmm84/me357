%% Example 1
clc; clear; close all;

H = tf([1 3 4],[1 2 1]);    % make TF
step(H);
grid on;

%% Example 2
clc; clear; close all;

H = tf([3 4],[1 2 4]);      % make TF
step(H);
grid on;

%% Example 3
clc; clear; close all;

H = tf([1 1],[1 7 16 12]);  % make TF
step(H);
grid on;

%% Example 4
clc; clear; close all;

H = tf([1],[1 5 11 15]); %#ok<*NBRAK>
step(H);
grid on;