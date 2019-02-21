fs = 44100; % sample rate

lengthSound = 5; % in [s]
totNumSamples = lengthSound * fs;

x = zeros(totNumSamples, 1);
y = zeros(totNumSamples, 1);

LP = zeros(totNumSamples, 1);
AP = zeros(totNumSamples, 1);

excitationL = fs*2; % excitation length

x(1:excitationL) = rand(excitationL, 1) * 2 - 1;
freq = 441;
g = 0;

pluckPos = 0.5;
DL1Length = floor(0.5 * fs / freq);

DL1 = zeros(DL1Length, 1);
DL2 = zeros(DL1Length, 1);
freq = 438.4;
x(1:excitationL) = sawtooth(2*pi*freq*[0:1/fs:(excitationL-1)/fs]);
pickupPos = 0.7;
for n = 2:totNumSamples
    
    %add the input
    DL1(floor(DLLength * pluckPos)) = DL1(floor(DLLength * pluckPos)) + x(n);
    DL2(floor(DLLength * pluckPos)) = DL2(floor(DLLength * pluckPos)) + x(n);
    
    % lowpass
    LP(n) = (DL2(1) + DL2(2)) / 2;
    
    % allpass
    AP(n) = -g * LP(n) + LP(n - 1) + g * AP(n - 1);
    
    %move samples around
    DL2(1:end - 1) = DL2(2:end);
    DL2(end) = -DL1(end);
    DL1(2:end) = DL1(1:end-1);
    DL1(1) = -AP(n);
    
    %obtain output
    y(n) = DL1(pickupPos * DLLength) + DL2(pickupPos * DLLength); 
    
end

%plot output
plot(y);