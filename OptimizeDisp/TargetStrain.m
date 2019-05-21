function error = TargetStrain(a, Rendo, Repi, Zbot, Ztop, thetaEndo, thetaMid, thetaEpi, TargetStrains, Weights, Phi, betamax, OptBeta, DebugFlag)
                 
% Deformation is axialsymmetric so the value of phi does not really matter in this case
error = 0.0;

NsamplePointsR = 2;
NsamplePointsZ = 3;

R = linspace(Rendo, Repi, NsamplePointsR);
Z = linspace(Zbot,  Ztop, NsamplePointsZ);

for i = 1:NsamplePointsR
    for j = 1:NsamplePointsZ
        
        [xa, ya, za, Fan] = AnalyticalF(Phi, R(i), Z(j), Zbot, Ztop, a, betamax, OptBeta);
        
        Can = transpose(Fan)*Fan;
        
        eR = [ cos(Phi); sin(Phi); 0.0];
        eT = [-sin(Phi); cos(Phi); 0.0];
        eZ = [ 0.0; 0.0; 1.0];
        fiber = Microstructure(thetaEndo, thetaMid, thetaEpi, Rendo, Repi, R(i), eR, eT, DebugFlag);

        J = det(Fan);
        ELL = 0.5*(eZ'*Can*eZ - 1.0);
        ECC = 0.5*(eT'*Can*eT - 1.0);
        ERR = 0.5*(eR'*Can*eR - 1.0);
        EFF = 0.5*(fiber(1:3)*Can*transpose(fiber(1:3)) - 1.0);
        
        error = error + Weights(1)  *(J   - TargetStrains(1)  )^2 + .... % J   term
                        Weights(2)  *(ELL - TargetStrains(2)  )^2 + .... % ELL term
                        Weights(2+i)*(ECC - TargetStrains(2+i))^2 + .... % ECC term
                        Weights(4+i)*(ERR - TargetStrains(4+i))^2 + .... % ERR term
                        Weights(7)  *(EFF - TargetStrains(7)  )^2;       % EFF term
    end
end











