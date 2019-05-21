function [Invariants, InvariantNames] = ComputeInvariants(F, Fiber, Lvec, Cvec, Rvec)
                                        
fib = transpose(Fiber(1:3));

C = transpose(F)*F; % right Cauchy Green deformation tensor
E = 0.5*(C - eye(3,3));   % Green-Lagrangian strain tensor
% M = fib*transpose(fib); 

% Compute invariants of the right Cauchy Green deformation tensor
Invariants(1) = trace(C);                           InvariantNames(1,:) = 'I1_'; % I1 tr(C)
Invariants(2) = 0.5*( trace(C)^2 - trace(C*C));     InvariantNames(2,:) = 'I2_'; % I2 1/2 (tr(C)^2 - tr(C*C))
Invariants(3) = det(C);                             InvariantNames(3,:) = 'I3_'; % I3 det(C)
Invariants(4) = transpose(fib)*C*fib;               InvariantNames(4,:) = 'I4_'; % I4
Invariants(5) = transpose(fib)*E*fib;               InvariantNames(5,:) = 'EFF'; % Fiber strain
Invariants(6) = transpose(Lvec)*E*Lvec;             InvariantNames(6,:) = 'ELL'; % longitundinal strain
Invariants(7) = transpose(Cvec)*E*Cvec;             InvariantNames(7,:) = 'ECC'; % circumferential strain
Invariants(8) = transpose(Rvec)*E*Rvec;             InvariantNames(8,:) = 'ERR'; % radial strain

end












