function afinal = OptimizeDisp(Rendo, Repi, Zbot, Ztop, thetaEndo, thetaMid, thetaEpi, TargetStrains, Weights, betamax, OptBeta, DebugFlag)

Phi    = 0.0; % It decides which point around circumference to probe - It does not matter since the phantom is axial-symmetric



%% Optimize displacement parameters to obtain correct cardiac deformation
switch OptBeta
    case 0 % Assign a prescribed twist
        a0 = [0.0, 0.0, 0.0, 0.0];
    case 1 % Compute beta as well
        a0 = [0.0, 0.0, 0.0, 0.0, 0.0];  
end

options = optimset('MaxFunEvals',10^5, 'MaxIter',10^5, 'PlotFcns', @optimplotfval, 'TolFun', 10^-7, 'TolX', 10^-7);
[afinal, fval, exitflag] = fminsearch(@(a) TargetStrain(a, Rendo, Repi, Zbot, Ztop, ....
                                                        thetaEndo, thetaMid, thetaEpi, ....
                                                        TargetStrains, Weights, Phi, betamax, OptBeta, DebugFlag), a0, options);

if (DebugFlag == 1)
    disp('Optimization results');
    disp(['fval     = ', num2str(fval)]);
    disp(['exitflag = ', num2str(exitflag)]);
    disp('                 J        ELL       ECC       ERR       EFF');

    %% Analyze strain values
    NsamplePointsR = 2;
    NsamplePointsZ = 3;
    Location = ['Endo - Zbot'; ....
                'Endo - Zmid'; ....
                'Endo - Ztop'; ....
                'Epi  - Zbot'; ....
                'Epi  - Zmid'; ....
                'Epi  - Ztop'];
    R = linspace(Rendo, Repi, NsamplePointsR);
    Z = linspace(Zbot,  Ztop, NsamplePointsZ);

    % Location = ['Endo', 'Zbot'; 'Epi ', 'Ztop'];
    for i = 1:NsamplePointsR
        for j = 1:NsamplePointsZ
            [xa, ya, za, Fan] = AnalyticalF(Phi, R(i), Z(j), Zbot, Ztop, afinal, betamax, OptBeta);

            eR = [ cos(Phi); sin(Phi); 0.0];
            eT = [-sin(Phi); cos(Phi); 0.0];
            eZ = [ 0.0; 0.0; 1.0];

            C = transpose(Fan)*Fan;

            J = det(Fan);
            ELL = 0.5*(eZ'*C*eZ - 1.0);
            ECC = 0.5*(eT'*C*eT - 1.0);
            ERR = 0.5*(eR'*C*eR - 1.0);
            fiber = Microstructure(thetaEndo, thetaMid, thetaEpi, Rendo, Repi, R(i), eR, eT, DebugFlag);
            EFF = 0.5*(fiber(1:3)*C*transpose(fiber(1:3)) - 1.0);
            
            disp([Location((i-1)*NsamplePointsZ + j, :),'   ',num2str(J),'   ',num2str(ELL),'   ',num2str(ECC),'   ',num2str(ERR),'   ',num2str(EFF)]);
        end
    end
    
    l_over_L = 1+afinal(4);
    rendo = afinal(1) + Rendo*(1+afinal(2)) + afinal(3)*Rendo^2;
    
    EF = 1.0 - l_over_L*(rendo/Rendo)^2;
    disp('');
    disp(['Ejection fraction = ', num2str(EF)]);
    disp('');
    afinal
end



end % function


