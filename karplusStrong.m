fs = 44100; % sample rate

lengthSound = 5; % in [s]
totNumSamples = lengthSound * fs;

x = zeros(totNumSamples, 1);
y = zeros(totNumSamples, 1);
LP = zeros(totNumSamples, 1);
AP = zeros(totNumSamples, 1);
x(1:100) = rand(100, 1) * 2 - 1;
freq = 440;
g = 0.1;
N = floor(fs / freq);
for n = N+2:totNumSamples
    LP(n) = (y(n - N) + y(n - N - 1)) / 2;
    AP(n) = -g * LP(n) + LP(n - 1) + g * AP(n - 1);
    y(n) = x(n - N) + AP(n);
end

plot(y);