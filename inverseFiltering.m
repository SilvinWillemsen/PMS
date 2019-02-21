fs = 44100;
bw = 5;
R = exp( - pi * bw / fs);
xylo = audioread("xylophone.wav");
output = xylo;
for mode = 1:10
    test = abs(fft(output));
    freq = find(test(1:floor(length(output)/2)) == max(test(1:floor(length(output)/2))));
    freq = (freq / length(output)) * fs;
    z = R * exp(j * 2 * pi * freq / fs);
    B = [1, -(z + conj(z)), z * conj(z)];
    r = 0.9;
    A = B .* (r.^[0 : length(B) - 1]);
    output = filter(B,A, output);
end
nonZero = find(output ~= 0);
nonZero = nonZero(1);
audiowrite("excitation.wav", output(nonZero:nonZero+100), fs) 
plot(output);