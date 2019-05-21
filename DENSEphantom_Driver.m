% % Luigi Perotti and Ilya Verzhbinsky - 2018 %
% DENSE PHANTOM DRIVER CODE
% 
% This script computes an analytical cardaic-like displacement field and 
% simulates 3D DENSE imaging data.

clear ;
close all ;
clc;
addpath(genpath('./'));

StartTime = cputime;

%% Input parameters

% Output path for vtk images and displacement files
OutputPath = ''; 

OutputPathVTK = fullfile(OutputPath,'VTK');      %Export path for vtk files
OutputPathMat = fullfile(OutputPath,'Matfiles'); %Export path for mat files
OutputPathIm  = fullfile(OutputPath,'Images');   %Export path for images

if ~isdir(OutputPathVTK); mkdir(OutputPathVTK); end
if ~isdir(OutputPathMat); mkdir(OutputPathMat); end
if ~isdir(OutputPathIm);  mkdir(OutputPathIm);  end

DebugFlag = 1;  % Print extra information



% Geometry of analytical (cylindrical) phantom
Rendo  = 25.0; % [mm]
Repi   = 35.0; % [mm]
Rmin   = 5.0;  % [mm] % Displacement field is singular at (0,0) - Rmin conservatively smaller than expected Rendo 
Zbot   = 0.0;  % [mm]
Ztop   = 16.0; % [mm]



% Model based myofiber orientation - Needed if the analytical phantom parameters are computed to obtain a certain myofiber stretch
thetaEndo =  37.0*pi/180.0;     % Endo fiber angle in E1     - (Journal of Biomechanics 41 (2008) 3219-3224) % 37 in E2  % 52 in E1
thetaMid  = -9.0*pi/180.0;      % Mid wall fiber angle in E1 - (Journal of Biomechanics 41 (2008) 3219-3224) % -9 in E2  % -6 in E1
thetaEpi  = -45*pi/180.0;       % Epi fiber angle in E1      - (Journal of Biomechanics 41 (2008) 3219-3224) % -45 in E2 % -49 in E1



% Define cardiac motion
ComputePhantom = 1; % 0 to use predefined values of a, 1 to recompute parameters of analytical phantom
afinal = [-23.0423    0.9233   -0.0099   -0.1145    1.5361]; % Used if ComputePhantom = 0
TargetStrains = [1.0,  -0.15, -0.2, -0.16, 0.45, 0.30, -0.15 ]; % J ELL ECC_endo ECC_epi ERR_endo ERR_epi EFF
Weights       = [0.5,   0.2 ,  0.5,   0.5,  0.1,  0.1,  1.0  ];



OptBeta   = 1;      % 0 -> fixed torsion; 1 -> minimize parameters to compute torsion
betamax = 0*pi/180; % [rad] Max torsion - Used if OptBeta = 0



TimeSteps = 40;
alphaTime = sin([0:pi/(TimeSteps-1):pi]); % Time variation of displacement magnitude. 
% It can be anything and different parts of the displacement field can also have their own time variation



% Define sequence parameters
% Field o view
Xlim = [-40, 40];    % [mm]
Ylim = [-40, 40];    % [mm]

zconst  = 4;        % Choose slice location in Z % [mm]
hx      = 2.5;      % In plane X grid size % [mm]
hy      = 2.5;      % In plane Y grid size % [mm]
hz      = 8.0;      % Through plane voxel size % [mm]
sampleX = 2;        % Sample points per voxel in X
sampleY = 2;        % Sample points per voxel in Y
sampleZ = 3;        % Sample points per voxel in Z
ke_x    = 0.08;     % cycles/mm
ke_y    = 0.08;     % cycles/mm
ke_z    = 0.08;     % cycles/mm

% Noise parameters 
SNR = [2,5,10,20,30,40,60,80,160,320];     % Normalized magnitude is 1 (||[Mx,My,Mz|| = 1)
reps = 1; %Number of reptitions per SNR 


% Visualization only parameters
NpointZ = 6;   % Number of elements in Z direction used for vtk representation
NpointR = 3;   % Number of elements in radial direction used for vtk representation
NpointC = 65;  % Number of elements in circumferential direction used for vtk representation
PhaseElType = 21; % Linear quadrilateral elements
CellScalarsName = ['Magnitude';'Phase_X  ';'Phase_Y  ';'Phase_Z  '];
PhaseFileName = [OutputPathVTK,'PhaseData'];
PhantomMeshSize = 1.0; % Not used unless NpointZ, NpointR, NpointC are less than 1


%% Compute motion coefficients for analytical phantom
if (ComputePhantom == 1)
    afinal = OptimizeDisp(Rendo, Repi, Zbot, Ztop, ....
                          thetaEndo, thetaMid, thetaEpi, ....
                          TargetStrains, Weights, ....
                          betamax, OptBeta, DebugFlag);
end

                  
                  

%% Visualize phan tom displacement in vtk
AnalyticalPhantom(afinal, betamax, OptBeta, Rendo, Repi, Zbot, Ztop, PhantomMeshSize, ....
                  NpointZ, NpointR, NpointC, ....
                  thetaEndo, thetaMid, thetaEpi, ....
                  alphaTime, OutputPathVTK, DebugFlag);
                
close                 
                  
                  

%% Generate phase data from displacement data
           
% Define solver options to transform Lagrangian to Eulerian displacements
tol = 5.0e-4; % This should be lower - Using the Optimization Toolbox is it possible
options = optimset('MaxFunEvals',10000,'MaxIter',10000);

% Generate grid points
gridX = [Xlim(1): hx :Xlim(2)]; NvoxelX = length(gridX)-1;
gridY = [Ylim(1): hy :Ylim(2)]; NvoxelY = length(gridY)-1;
[xPhase, connPhase] = PhaseMesh(gridX, gridY, zconst); % Element nodal coordinates and connectivity to plot phase data in vtk

Zg = linspace(zconst - 0.5*hz, zconst + 0.5*hz, sampleZ*2+1);
Zg(1:2:length(Zg)) = [];

NpointsPerVoxel = sampleX*sampleY*sampleZ;

% Loop over time (1 image per time)
for s = 20%1:length(alphaTime) % loop over time
    % Initialize Disp array to zero
    dispAllX = zeros(NvoxelX*sampleX, NvoxelY*sampleY, sampleZ); % Only needed if this information need to be saved for postprocessing
    dispAllY = zeros(NvoxelX*sampleX, NvoxelY*sampleY, sampleZ); % Only needed if this information need to be saved for postprocessing
    dispAllZ = zeros(NvoxelX*sampleX, NvoxelY*sampleY, sampleZ); % Only needed if this information need to be saved for postprocessing
    
    % Initialize phase and magnitude arrays to zero
    magnitude = zeros(NvoxelX, NvoxelY);
    phaseX    = zeros(NvoxelX, NvoxelY);
    phaseY    = zeros(NvoxelX, NvoxelY);
    phaseZ    = zeros(NvoxelX, NvoxelY);
    
    for i = 1:NvoxelX        % Loop over voxels
        Xg = linspace(gridX(i), gridX(i+1), sampleX*2+1);
        Xg(1:2:length(Xg)) = [];
        for j = 1:NvoxelY
            Yg = linspace(gridY(j), gridY(j+1), sampleY*2+1);
            Yg(1:2:length(Yg)) = [];
            % Loop over multiple points in the voxel to consider intravoxel dephasing
            NphaseToAverage = 0;
            for is = 1:sampleX
                for js = 1:sampleY
                    for ks = 1:sampleZ
                        xgiven = [Xg(is), Yg(js), Zg(ks)];  % Location in Eulerian coordinates for which we want to compute the displacement
                        if ( norm(xgiven(1:2)) > Rmin )
                            indX = (i-1)*sampleX + is;
                            indY = (j-1)*sampleY + js;
                            indZ = ks;
                            X0 = xgiven - [dispAllX(indX,indY,indZ), dispAllY(indX,indY,indZ), dispAllZ(indX,indY,indZ)];    % Initial guess
                            
                            a_s = afinal * alphaTime(s);
                            beta_s = betamax * alphaTime(s);
                            
                            [Xval, fval] = fsolve(@(X) EulerianMap(X, xgiven, Zbot, Ztop, a_s, beta_s, OptBeta), X0, options); % Find X corresponding to a given x
                            
                            % If you don't have the optimization toolbox,
                            % use the following:
                            % [Xval, fval] = fminsearch(@(X) EulerianMapNoOptToolbox(X, xgiven, Zbot, Ztop, a_s, beta_s, OptBeta), X0, options); % Find X corresponding to a given x
                            
                            
                            if ( norm(fval) < tol && IsMyocardium(Xval, Rendo, Repi, Zbot, Ztop) == 1)
                                dispXsingle = xgiven(1) - Xval(1);
                                dispYsingle = xgiven(2) - Xval(2);
                                dispZsingle = xgiven(3) - Xval(3);
                                
                                dispAllX(indX,indY,indZ) = dispXsingle;
                                dispAllY(indX,indY,indZ) = dispYsingle;
                                dispAllZ(indX,indY,indZ) = dispZsingle;
                                
                                phaseX(i,j) = phaseX(i,j) + dispXsingle;
                                phaseY(i,j) = phaseY(i,j) + dispYsingle;
                                phaseZ(i,j) = phaseZ(i,j) + dispZsingle;
                                NphaseToAverage = NphaseToAverage+1;
                                
                                magnitude(i,j) = magnitude(i,j) + 1.0;
                            else
                                dispAllX(indX,indY,indZ) = 0;
                                dispAllY(indX,indY,indZ) = 0;
                                dispAllZ(indX,indY,indZ) = 0;
                            end
                        end 
                    end % loop over ks
                end % loop over js
            end % loop over is
            
            % Average displacement over voxel - Compute phase and magnitude
            % Compute magnitude 
            magnitude(i,j) = magnitude(i,j)/NpointsPerVoxel;
            
            if (NphaseToAverage == 0)
                NphaseToAverage = 1;
            end
            
            phaseX(i,j) = ( phaseX(i,j) * ke_x / NphaseToAverage);
            phaseY(i,j) = ( phaseY(i,j) * ke_y / NphaseToAverage);
            phaseZ(i,j) = ( phaseZ(i,j) * ke_z / NphaseToAverage);
            
        end % loop over j
    end % loop over i
    
    for rep = 1:reps
        for n = SNR % loop through noise
            
            % + 0.5 so that it will wrap after +- pi ; - 0.5 so that zero displacement corresponds to zero phase
            phaseX = mod(phaseX + 0.5, 1) - 0.5;
            phaseY = mod(phaseY + 0.5, 1) - 0.5; 
            phaseZ = mod(phaseZ + 0.5, 1) - 0.5;

            % Add noise to magnitude and phase data
            [magnitude_wN, phaseX_wN, phaseY_wN, phaseZ_wN] = AddNoiseToData_DENSE(magnitude, phaseX, phaseY, phaseZ, n);


            phaseX_wN = mod(phaseX_wN + 0.5, 1) - 0.5;
            phaseY_wN = mod(phaseY_wN + 0.5, 1) - 0.5;
            phaseZ_wN = mod(phaseZ_wN + 0.5, 1) - 0.5;

            CellScalars(:,1) = reshape(magnitude_wN, NvoxelX*NvoxelY, 1);
            CellScalars(:,2) = reshape(phaseX_wN, NvoxelX*NvoxelY, 1);
            CellScalars(:,3) = reshape(phaseY_wN, NvoxelX*NvoxelY, 1);
            CellScalars(:,4) = reshape(phaseZ_wN, NvoxelX*NvoxelY, 1);



            % Plot phase data to VTK
            PlotToVTK(xPhase, connPhase-1, PhaseElType, [], [], CellScalars, CellScalarsName, PhaseFileName, s);



            % Plot images directly in matlab - debug only
            figure(10);
            subplot(2,2,1)
            pcolor(phaseX_wN);
            title(['X [cycles]  ', num2str(s)]); xlabel('X'); ylabel('Y'); caxis([-0.5 0.5]); colormap gray; colorbar; axis square; set(gca,'fontsize',16);
            
            subplot(2,2,2)
            pcolor(phaseY_wN);
            title(['Y [cycles]  ', num2str(s)]); xlabel('X'); ylabel('Y'); caxis([-0.5 0.5]); colormap gray; colorbar; axis square; set(gca,'fontsize',16);
            
            subplot(2,2,3)
            pcolor(phaseZ_wN);
            title(['Z [cycles]  ', num2str(s)]); xlabel('X'); ylabel('Y'); caxis([-0.5 0.5]); colormap gray; colorbar; axis square; set(gca,'fontsize',16);
            
            subplot(2,2,4)
            pcolor(magnitude_wN);
            title(['Magnitude  ', num2str(s)]); xlabel('X'); ylabel('Y'); caxis([0 1]); colormap gray; colorbar; axis square; set(gca,'fontsize',16);
            pause(0.1)

            set(gcf, 'Position', [100, 100, 1000, 1000])

            % Save DispX DispY DispZ to file
            OutputPathMat_temp = fullfile(OutputPathMat,['Zpos_',num2str(zconst)],['SNR_',num2str(n,'%02d')],['Rep_',num2str(rep)]);
            if ~isdir(OutputPathMat_temp); mkdir(OutputPathMat_temp); end

            save(fullfile(OutputPathMat_temp,['phaseX_',num2str(s,'%02d'),'.mat']),'phaseX_wN');
            save(fullfile(OutputPathMat_temp,['phaseY_',num2str(s,'%02d'),'.mat']),'phaseY_wN');
            save(fullfile(OutputPathMat_temp,['phaseZ_',num2str(s,'%02d'),'.mat']),'phaseZ_wN');
            save(fullfile(OutputPathMat_temp,['mag_',num2str(s,'%02d'),'.mat']),'magnitude_wN');

        end
    end

end

disp(['Elapsed time = ',num2str(cputime - StartTime)]);
    	
%%
save(fullfile(OutputPath,'PhantomData.mat'))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    









