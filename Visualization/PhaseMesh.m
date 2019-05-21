function [xPhase, connPhase] = PhaseMesh(gridX, gridY, zconst)

ind = 1;
Nx = length(gridX);
Ny = length(gridY);
for i = 1:Nx
    for j = 1:Ny
        xPhase(ind, 1:3) = [gridX(i), gridY(j), zconst];
        ind = ind +1;
    end
end

ind = 1;
for i = 1:Nx-1
    for j = 1:Ny-1
        connPhase(ind, 1:4) = [(i-1)*Nx+j , i*Nx+j, i*Nx+j+1, (i-1)*Nx+j+1];
        ind = ind +1;
    end
end









