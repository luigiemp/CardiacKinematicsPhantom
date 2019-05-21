function [nodesCart, nodesPolar, conn] = PhantomHexaMesh(Rendo, Repi, Zbot, Ztop, h, NpointZ, NpointR, NpointC)

% Create cylindrical mesh based on the (Rendo, Repi, Zbot, Ztop, h)
% Return nodes in cartesian and polar coordinates together with
% connectivity table



%% Generate nodal coordinates
DeltaH = Ztop - Zbot;

if (NpointZ < 1)
    NpointZ = ceil(DeltaH/h);
end

if (NpointR < 1)
    NpointR = ceil((Repi - Rendo)/h);
end
if (NpointC < 1)
    NpointC = ceil(pi*(Rendo + Repi)/h);
end

ind = 1;
for k = 0:NpointZ;
    z = Zbot + k*DeltaH/NpointZ;
    for j = 0:NpointR;
        r = Rendo + (Repi-Rendo)*j/NpointR;
        for i = 1:NpointC ;
            theta = 2*pi*i/NpointC;
            nodesCart(ind,  1:3) = [r*cos(theta), r*sin(theta), z];
            nodesPolar(ind, 1:3) = [theta, r, z];
            ind = ind + 1;
        end
    end
end



%% Generate slice connectivity - Formatted for VTK!!!
ind = 1;

for k = 0:NpointZ-1
    upZ = NpointC*(NpointR+1);
    for j = 0:NpointR-1
        upR = NpointC;
        for i = 1:NpointC-1 
            conn(ind, 1:8) = [i, i+upR, i+1+upR, i+1, ....
                              i+upZ, i+upR+upZ, i+1+upR+upZ, i+1+upZ] + upR*j + upZ*k;
            ind = ind + 1;
        end
        conn(ind, 1:8) = [NpointC, NpointC+upR, 1+upR, 1, ....
                          NpointC+upZ, NpointC+upR+upZ, 1+upR+upZ, 1+upZ] + upR*j + upZ*k;
        ind = ind + 1;
    end
    
end









