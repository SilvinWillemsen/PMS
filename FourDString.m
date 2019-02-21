clear all;
close all;

fs = 44100;
k = 1/fs;
numMasses = 50;
m = 0.01;
K = 1000000;
K3 = 000000;
z = 0.0;
numDimensions = 4;
uNext = zeros(numMasses, numDimensions);
u = zeros(numMasses, numDimensions);
width = 14;
posY = floor(1/3 * numMasses);
posZ = floor(1/5 * numMasses);
pos4 = floor (4/7 * numMasses);
% u(posY - width / 2 : posY + width / 2, 2) = sin(pi * [0:width] / width);
% u(posZ - width / 2 : posZ + width / 2, 3) = sin(2 * pi * [0:width] / width);
% u(:,1) = 1/numMasses:1/numMasses:1;
% u(pos4, 4) = 1;
l0 = 1/(10 * numMasses); % resting length of the springs
u(:,1) = 1/numMasses : 1/numMasses : 1;
% u(:,4) = 1/numMasses : 1/numMasses : 1;
u(floor(numMasses / 2),2) = 1;
uPrev = u;
uNext = u;

% figure;
ylim([-1 1]);
F1 = zeros(1, numDimensions);
F2 = zeros(1, numDimensions);
M = m / k^2;
Z = z / k;
sigma = 0 / k;
lengthSound = fs;
output = zeros(lengthSound, 1);

zeroVec = zeros(numMasses, 1);
blueVec = ones(numMasses, 1);
outIdx = floor(numMasses / 4);
for n = 1:lengthSound
%     for i = 2:numMasses-1
%         
%         dist = sqrt(sum((u(i,:) - u(i-1,:)).^2));
%         prevDist= sqrt(sum((uPrev(i,:) - uPrev(i-1,:)).^2));
%         Ftot = -K * (dist - l0) - K3 * (dist - l0)^3 - Z * (dist - prevDist);
%         
%         dist2 = sqrt(sum((u(i,:) - u(i+1,:)).^2));
%         prevDist2 = sqrt(sum((uPrev(i,:) - uPrev(i+1,:)).^2));
%         Ftot2 = -K * (dist2 - l0) - K3 * (dist2 - l0)^3 - Z * (dist2 - prevDist2);
%         
%         F = Ftot * (u(i,:) - u(i-1,:)) / dist;
%         F2 = Ftot2 * (u(i,:) - u(i+1,:)) / dist2;
% %         if i ~= numMasses
%             uNext(i,:) = 2 * u(i,:) - uPrev(i,:) + F / M;
%             uNext(i,:) = uNext(i,:) + F2 / M;
% %         end
% %         if i ~= 2
% %             uNext(i-1, :) = uNext(i-1, :) + F / M;
% %         end
%     end
    for i = 2:numMasses
        
        dist = sqrt(sum((u(i,:) - u(i-1,:)).^2));
        prevDist= sqrt(sum((uPrev(i,:) - uPrev(i-1,:)).^2));
        Ftot = -K * (dist - l0) - K3 * (dist - l0)^3 - Z * (dist - prevDist);

        F = Ftot * (u(i,:) - u(i-1,:)) / dist;
        if i ~= numMasses
            uNext(i,:) = (2 - sigma / M) * u(i,:) + (-1 + sigma / M) * uPrev(i,:) + F / M;
        end
        if i ~= 2
            uNext(i-1, :) = uNext(i-1, :) - F / M;
        end
    end
    output(n,1) = u(outIdx, 1) - outIdx/numMasses;
    output(n,2) = u(outIdx, 2);
    output(n,3) = u(outIdx, 3);
    output(n,4) = u(outIdx, 4) - outIdx/numMasses;
%     output(n,5) = sum([x, u(outIdx,2:end)]);

    uPrev = u;
    u = uNext;
    
    if mod(n, 1) == 0 %&& n > lengthSound
        subplot(3,1,1)
        scatter3(u(:,1) * numMasses, u(:,2), u(:,3), 400, [u(:,4) + 0.5 - [1/numMasses : 1/numMasses : 1]', zeroVec, blueVec], '.');
        xlim([0 numMasses]);
        ylim([-1 1])
        zlim([-1 1])
        subplot(3,1,2)
        plot(u(:,4)- [1/numMasses : 1/numMasses : 1]')
        subplot(3,1,3)
        plot(u(:,1)- [1/numMasses : 1/numMasses : 1]')
        drawnow;
    end
end
plot(output(:,2))