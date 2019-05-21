function [xa, ya, za, Fan] = AnalyticalF(Phi, R, Z, Zbot, Ztop, a, beta, OptBeta)

% Return in plane deformation gradient F
vR = a(1) + R*(1+a(2)) + a(3)*R^2;
vZ = Z + a(4)*Z;
switch OptBeta
    case 0
        vT = 0.0;
        dvT_dZ = 0.0;
    case 1
        vT = a(5)*(Z-Zbot)/(Ztop-Zbot);
        dvT_dZ = a(5)/(Ztop-Zbot);
end

dvR_dR = 1 + a(2) + 2*a(3)*R;
dvR_dT = 0.0;
dvR_dZ = 0.0;

dvT_dR = 0.0;
dvT_dT = 0.0;

dvZ_dR = 0.0;
dvZ_dT = 0.0;
dvZ_dZ = 1.0 + a(4);

Fp = [dvR_dR,     (dvR_dT - vT)/R,   dvR_dZ; ....
      dvT_dR,     (dvT_dT + vR)/R,   dvT_dZ; ....
      dvZ_dR,     (dvZ_dT)/R,        dvZ_dZ];


eR = [ cos(Phi); sin(Phi); 0.0];
eT = [-sin(Phi); cos(Phi); 0.0];
eZ = [ 0.0; 0.0; 1.0];


F_inplane = Fp(1,1) * eR*eR' + Fp(1,2) * eR*eT' + Fp(1,3) * eR*eZ' + ....
            Fp(2,1) * eT*eR' + Fp(2,2) * eT*eT' + Fp(2,3) * eT*eZ' + ....
            Fp(3,1) * eZ*eR' + Fp(3,2) * eZ*eT' + Fp(3,3) * eZ*eZ';

switch OptBeta
    case 0
        beta = beta*(Z-Zbot)/(Ztop-Zbot);
        dbetadR = 0;
        dbetadZ = beta/(Ztop-Zbot);

        Frot = [    cos(beta)-R*sin(Phi+beta)*cos(Phi)*dbetadR, -sin(beta)-R*sin(Phi+beta)*sin(Phi)*dbetadR, -R*sin(Phi+beta)*dbetadZ; ....
                    sin(beta)+R*cos(Phi+beta)*cos(Phi)*dbetadR,  cos(beta)+R*cos(Phi+beta)*sin(Phi)*dbetadR,  R*cos(Phi+beta)*dbetadZ; ....
                    0.0,                                        0.0,                                          1.0];

        Fan = F_inplane * Frot;
        
    case 1
        Fan = F_inplane;
end

% Transform vR, vT, and vZ from polar to cartesian coordinates
xa = vR*cos(Phi) - vT*sin(Phi);
ya = vR*sin(Phi) + vT*cos(Phi);
za = vZ;

end

