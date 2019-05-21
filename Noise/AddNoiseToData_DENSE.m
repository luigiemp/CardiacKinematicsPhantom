
% Simulate DENSE x-y-z phase assuming a balanced 4 point encoding scheme.

% Zhong, X., Helm, P. A., & Epstein, F. H. (2009). Balanced multipoint 
% displacement encoding for DENSE MRI. Magnetic Resonance in Medicine: 
% An Official Journal of the International Society for Magnetic Resonance 
% in Medicine, 61(4), 981-988.



function [magnitude_wN, phaseX_wN, phaseY_wN, phaseZ_wN] = AddNoiseToData_DENSE(magnitude, phaseX, phaseY, phaseZ, SNR)

[Nx, Ny] = size(magnitude);

sqrt3over3 = sqrt(3.0)/3.0;
A = [-sqrt3over3, -sqrt3over3, -sqrt3over3, 1; ....
      sqrt3over3,  sqrt3over3, -sqrt3over3, 1; ....
      sqrt3over3, -sqrt3over3,  sqrt3over3, 1; ....
     -sqrt3over3,  sqrt3over3,  sqrt3over3, 1];
 
sqrt3over4 = sqrt(3.0)/4.0;
Ainv = [-sqrt3over4,  sqrt3over4,  sqrt3over4, -sqrt3over4; ....
        -sqrt3over4,  sqrt3over4, -sqrt3over4,  sqrt3over4; ....
        -sqrt3over4, -sqrt3over4,  sqrt3over4,  sqrt3over4; ....
         0.25,  0.25,  0.25,  0.25];
     
NoiseStandardDeviation = 2.0/SNR; % Because we are averaging 4 normally distributed noise and we want 
% the summed noise to have standard deviation equal to 1/SNR

S0_noise = NoiseStandardDeviation*complex(randn(Nx, Ny), randn(Nx, Ny));
S1_noise = NoiseStandardDeviation*complex(randn(Nx, Ny), randn(Nx, Ny));
S2_noise = NoiseStandardDeviation*complex(randn(Nx, Ny), randn(Nx, Ny));
S3_noise = NoiseStandardDeviation*complex(randn(Nx, Ny), randn(Nx, Ny));

magnitude_wN = zeros(Nx, Ny);
phaseX_wN = zeros(Nx, Ny);
phaseY_wN = zeros(Nx, Ny);
phaseZ_wN = zeros(Nx, Ny);

for a = 1:Nx
    for b = 1:Ny
        
        phi_ab = 2*pi*A*[phaseX(a,b); phaseY(a,b); phaseZ(a,b); 0.0]; % Use background phase phi_b = 0.0
        
        S_ab = magnitude(a,b)*complex(cos(phi_ab), sin(phi_ab)); 
        
        % Compute phase with noise
        Swn = S_ab + [S0_noise(a,b); S1_noise(a,b); S2_noise(a,b); S3_noise(a,b)];
        
        phi_ab_wn = angle( Swn );
        
        phi_ab_wn = round( (phi_ab - phi_ab_wn)/(2*pi) ) + phi_ab_wn/(2*pi); % Unwrapping
       
        D_phi_b = Ainv*phi_ab_wn;
        % phi_b = D_phi_b(4);
            
        phaseX_wN(a,b) = D_phi_b(1); % - phi_b;
        phaseY_wN(a,b) = D_phi_b(2); % - phi_b;
        phaseZ_wN(a,b) = D_phi_b(3); % - phi_b;
        
        % Compute magnitude with noise
        magnitude_wN(a,b) = mean( abs( Swn ) );

    end
end










