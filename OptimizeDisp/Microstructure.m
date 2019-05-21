function fiber = Microstructure(ThetaEndo, ThetaMid, ThetaEpi, Rendo, Repi, R, eR, eT, DebugFlag)
    
    % Normalize alpha to [0, 1] endo to epi
    r = (R-Rendo)/(Repi-Rendo);
    
    % Compute angle to rotate fiber
    theta = 2.0 * ThetaEndo *(r - 0.5)*(r - 1.0) - ....
            4.0 * ThetaMid  *(r)*(r-1) + ....
            2.0 * ThetaEpi  *(r)*(r-0.5);
    
    % Compute fiber angle at element position 
    fiber(1:3) = eT*cos(theta) + cross(eR,eT)*sin(theta);
    
    % Check - This must be true by construction
    tol = 1.0e-8;
    if (DebugFlag == 1)
        if (dot(fiber(1:3), eR) > tol) % f0 and k0 are perpendicular
            dot(fiber(1:3), eR)
        end
        if (abs(norm(fiber(1:3)) - 1) > tol) % frot should be already normalized
            norm(fiber(1:3))
        end
    end
    fperp = cross(fiber(1:3), eR);
        
    fiber(4:6) = fperp;
    fiber(7:9) = eR;

end

