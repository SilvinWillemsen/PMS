% Backwards-time mass-spring-damping
% with multiple modes

clear all;
close all;
clc;

fs = 44100;
k = 1000;
f0 = 880;
m = 1;
B = 0.005;
numModes = 5;
excitation = audioread("excitation.wav");
lengthSound = fs;
excitation = [excitation; zeros(lengthSound - length(excitation), 1)];
 
for mode = 1:numModes
    f(mode) = mode * f0 * sqrt(1 + B * mode^2);
    k(mode) = m * (2 * pi * f(mode))^2;
    x(:, mode) = zeros(lengthSound,1);
    x(1:2, mode) = 0;
end

output = zeros(lengthSound,1);
R = 5;
T = 1/fs;

for n = 2:lengthSound - 1
%     x(n) = ((2 * m/T^2 + R / T) * x(n-1) - m/T^2 * x(n-2)) / (m/T^2 + R/T + k);
    for mode = 1 : numModes
        x(n+1, mode)= 2 * x(n, mode) - x(n - 1, mode) - (R/T * (x(n, mode) - x(n - 1, mode)) + k(mode) * x(n, mode)) / (m/T^2) + excitation(n); 
        output(n+1) = output(n+1) + x(n+1, mode);
    end
end
plot(output)

