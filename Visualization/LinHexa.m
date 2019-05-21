function [Fnum] = LinHexa(Xel, phiEl, numQP)

switch(numQP)
    case 0
        QP = [0.0, 0.0, 0.0];
    case 1
        qp = 1.0/sqrt(3);
        QP = [-qp, -qp, -qp; ....
               qp, -qp, -qp; ....
               qp,  qp, -qp; ....
              -qp,  qp, -qp; ....
              -qp, -qp,  qp; ....
               qp, -qp,  qp; ....
               qp,  qp,  qp; ....
              -qp,  qp,  qp];
end
  
for q=1:size(QP,1)
    % derivative of shape functions in isoparametric space (row 1 ->
    % derivative wrt x, row 2 -> derivative wrt y, row 3 derivative wrt z)
        DNr = 0.125*[(QP(q,2)-1)*(1-QP(q,3)), (1-QP(q,2))*(1-QP(q,3)), (1+QP(q,2))*(1-QP(q,3)), (1+QP(q,2))*(QP(q,3)-1), ....
                     (QP(q,2)-1)*(1+QP(q,3)), (1-QP(q,2))*(1+QP(q,3)), (1+QP(q,2))*(1+QP(q,3)),-(1+QP(q,2))*(1+QP(q,3)); ....
                     (QP(q,1)-1)*(1-QP(q,3)), (1+QP(q,1))*(QP(q,3)-1), (1+QP(q,1))*(1-QP(q,3)), (1-QP(q,1))*(1-QP(q,3)), ....
                     (QP(q,1)-1)*(1+QP(q,3)),-(QP(q,1)+1)*(1+QP(q,3)), (1+QP(q,1))*(1+QP(q,3)), (1-QP(q,1))*(1+QP(q,3)); ....
                     (QP(q,1)-1)*(1-QP(q,2)), (1+QP(q,1))*(QP(q,2)-1),-(1+QP(q,1))*(1+QP(q,2)), (QP(q,1)-1)*(1+QP(q,2)), ....
                     (1-QP(q,1))*(1-QP(q,2)), (1+QP(q,1))*(1-QP(q,2)), (1+QP(q,1))*(1+QP(q,2)), (1-QP(q,1))*(1+QP(q,2))];
        
        % Jacobian transformation matrix
        JT = DNr*Xel;
        detJ = det(JT);  % Element volume
        if (detJ < 0)
            disp('ERROR : element with negative Jacobian');
        end
        JTinv = inv(JT);
    
    DNx = JTinv*DNr; % Shape functions derivatives in  XYZ domain     
    
    
    % Compute deformation gradient F
    Fnum = zeros(3,3);
    for i = 1:8
        Fnum(1,1) = Fnum(1,1) + phiEl(i,1)*DNx(1,i);
        Fnum(1,2) = Fnum(1,2) + phiEl(i,1)*DNx(2,i);
        Fnum(1,3) = Fnum(1,3) + phiEl(i,1)*DNx(3,i);
        
        Fnum(2,1) = Fnum(2,1) + phiEl(i,2)*DNx(1,i);
        Fnum(2,2) = Fnum(2,2) + phiEl(i,2)*DNx(2,i);
        Fnum(2,3) = Fnum(2,3) + phiEl(i,2)*DNx(3,i);
        
        Fnum(3,1) = Fnum(3,1) + phiEl(i,3)*DNx(1,i);
        Fnum(3,2) = Fnum(3,2) + phiEl(i,3)*DNx(2,i);
        Fnum(3,3) = Fnum(3,3) + phiEl(i,3)*DNx(3,i);
    end % end F loop
    
end % end quadrature point loop


end % end Assemble Hexahedral element 
