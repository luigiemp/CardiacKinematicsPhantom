function [magnitude_wN, phaseX_wN, phaseY_wN, phaseZ_wN] = AddNoiseToData_Simple(magnitude, phaseX, phaseY, phaseZ, SNR)

[Nx, Ny] = size(magnitude);

NoiseStandardDeviation = sqrt(3)/SNR; % Because we are averaging 3 normally distributed noise and we want 
% the summed noise to have standard deviation equal to 1/SNR

Sx = magnitude.*complex(cos(phaseX*2*pi), sin(phaseX*2*pi));
Noise_x = NoiseStandardDeviation*(randn(Nx, Ny) + i*randn(Nx, Ny));
phaseX_wN = angle(Sx + Noise_x)/(2*pi); 
       
Sy = magnitude.*complex(cos(phaseY*2*pi), sin(phaseY*2*pi));
Noise_y = NoiseStandardDeviation*(randn(Nx, Ny) + i*randn(Nx, Ny));
phaseY_wN = angle(Sy + Noise_y)/(2*pi);

Sz = magnitude.*complex(cos(phaseZ*2*pi), sin(phaseZ*2*pi));
Noise_z = NoiseStandardDeviation*(randn(Nx, Ny) + i*randn(Nx, Ny));
phaseZ_wN = angle(Sz + Noise_z)/(2*pi);

magnitude_wN = (abs(Sx + Noise_x) + abs(Sy + Noise_y) + abs(Sz + Noise_z))/3.0;


% magnitude_Noise_X = NoiseStandardDeviation*randn(Nx, Ny);
% phaseNoiseX = -pi + 2*pi*rand(Nx, Ny);
% phaseX_wN = angle(magnitude.*complex(cos(phaseX*2*pi), sin(phaseX*2*pi)) + magnitude_Noise_X.*complex(cos(phaseNoiseX), sin(phaseNoiseX)))/(2*pi);
% 
% magnitude_Noise_Y = NoiseStandardDeviation*randn(Nx, Ny);
% phaseNoiseY = -pi + 2*pi*rand(Nx, Ny);
% phaseY_wN = angle(magnitude.*complex(cos(phaseY*2*pi), sin(phaseY*2*pi)) + magnitude_Noise_Y.*complex(cos(phaseNoiseY), sin(phaseNoiseY)))/(2*pi);
% 
% magnitude_Noise_Z = NoiseStandardDeviation*randn(Nx, Ny);
% phaseNoiseZ = -pi + 2*pi*rand(Nx, Ny);
% phaseZ_wN = angle(magnitude.*complex(cos(phaseZ*2*pi), sin(phaseZ*2*pi)) + magnitude_Noise_Z.*complex(cos(phaseNoiseZ), sin(phaseNoiseZ)))/(2*pi);









