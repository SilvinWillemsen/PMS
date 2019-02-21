clear all;
close all;

fs = 44100;
k = 1/fs;
Nx = 20;
Ny = 20;
m = 0.01;
K = 100000;
K3 = 000000;
z = 0.5;
sigma = 10 / k;
numDimensions = 3;
uNext = zeros(Nx * Ny, numDimensions);
u = zeros(Nx * Ny, numDimensions);
width = 14;
% posY = floor(1/3 * numMasses);
% posZ = floor(1/5 * numMasses);
% pos4 = floor (4/7 * numMasses);
% u(posY - width / 2 : posY + width / 2, 2) = sin(pi * [0:width] / width);
% u(posZ - width / 2 : posZ + width / 2, 3) = sin(2 * pi * [0:width] / width);
% u(:,1) = 1/numMasses:1/numMasses:1;
% u(pos4, 4) = 1;
l0 = 1/(8 * Nx); % resting length of the springs
for i = 1:Nx
    u([1:Ny] + (i - 1) * Ny, 1) = ones(Ny, 1) * i/Ny;
end
u1Init = u(:,1);
u(:, 2) = repmat([1/Ny : 1/Ny : 1]', Nx, 1);
% u(:,4) = 1/numMasses : 1/numMasses : 1;
u(Ny / 2 + Nx * Nx / 2, 3) = 1;
uPrev = u;
uNext = u;

% figure;
ylim([-1 1]);
F1 = zeros(1, numDimensions);
F2 = zeros(1, numDimensions);
M = m / k^2;
Z = z / k;
lengthSound = fs;
output = zeros(lengthSound, 1);

zeroVec = zeros(Nx * Ny, 1);
blueVec = ones(Nx * Ny, 1);
outIdx = floor(Ny / 2 + Nx * Nx / 2);

drawPlate = false;
for n = 1:lengthSound
%     vec = 2:numMasses;
    uNext = 2 * u - uPrev;
    
    for i = 2:Nx
        for j = 2:Ny
            idx = j + (i-1) * Nx;
            
            distX = sqrt(sum((u(idx, :) - u(idx - Ny, :)).^2));
            prevDistX = sqrt(sum((uPrev(idx, :) - uPrev(idx - Ny, :)).^2));
            FtotX = -K * (distX - l0) - K3 * (distX - l0)^3 - Z * (distX - prevDistX);
            
            distY = sqrt(sum((u(idx, :) - u(idx - 1, :)).^2));
            prevDistY = sqrt(sum((uPrev(idx, :) - uPrev(idx-1, :)).^2));
            FtotY = -K * (distY - l0) - K3 * (distY - l0)^3 - Z * (distY - prevDistY);

            FX = FtotX * (u(idx,:) - u(idx - Ny,:)) / distX;
            FY = FtotY * (u(idx,:) - u(idx-1,:)) / distY;
            if i ~= Nx
                if j ~= Ny
                    uNext(idx,:) = uNext(idx, :) + (FX + FY) / M;
                    if idx == length(u) - 5
                    end
                end
            end
            if j == 2 || i == Nx
            else
                uNext(idx - 1, :) = uNext(idx - 1, :) - FY / M;
            end
            if i == 2 || j == Ny
            else
                uNext(idx - Ny, :) = uNext(idx - Ny, :) - FX / M;
            end
        end
    end
    
    if drawPlate == true && mod(n, 10) == 0
        testMatX = reshape(u(:,1), Ny, Nx);
        testMatY = reshape(u(:,2), Ny, Nx);
        testMatZ = reshape(u(:,3), Ny, Nx);
        surf(testMatX, testMatY, testMatZ)
        view(60, 30)
        drawnow;
    end
    
    output(n,1) = u(outIdx, 1) - 1/2 - 1/Nx;
    output(n,2) = u(outIdx, 2) - 1/2 - 1/Ny;
    output(n,3) = u(outIdx, 3);
%     output(n,4) = u(outIdx, 4) - outIdx/Nx;
%     output(n,5) = sum([x, u(outIdx,2:end)]);

    uPrev = u;
    u = uNext;
end
plot(output(:,2))