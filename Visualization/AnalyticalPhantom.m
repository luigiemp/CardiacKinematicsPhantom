function AnalyticalPhantom(afinal, betamax, OptBeta, Rendo, Repi, Zbot, Ztop, h, NpointZ, NpointR, NpointC, thetaEndo, thetaMid, thetaEpi, alphaTime, OutputPathVTK, DebugFlag)

NumInv = 8;
Phi    = 0.0; % It decides which point to probe - It does not matter since the phantom is axial-symmetric

%% Compute initial nodal position, mesh, shape functions derivatives, and fibers
[nodesCart, nodesPolar, conn] = PhantomHexaMesh(Rendo, Repi, Zbot, Ztop, h, NpointZ, NpointR, NpointC);
NumDOF = size(nodesCart,1)*size(nodesCart,2);



%% Compute deformed heart along alphaTime

for s = 1:length(alphaTime)
    % Compute deformed heart at time s
    a_s = afinal * alphaTime(s);
    for i = 1:size(nodesPolar,1)
        [xa, ya, za, Fan] = AnalyticalF(nodesPolar(i,1), nodesPolar(i,2), nodesPolar(i,3), Zbot, Ztop, a_s, betamax*alphaTime(s), OptBeta);
        xt(i,1:3) = [xa, ya, za]; 
    end
    
    % Compute deformation gradient and plot its invariants - analytical phantom validation
	for e = 1:size(conn,1)
        
        % Using barycenter to compute analytical invariants
        Xel = nodesCart(conn(e,:),:);
        ElBar = mean(Xel);
        nodesQP(e,:) = mean(xt(conn(e,:),:));
        
        xqp = ElBar(1); yqp = ElBar(2); zqp = ElBar(3);
        rqp = sqrt(xqp^2 + yqp^2);
        
        Lvec = [0; 0; 1];
        Cvec = [yqp; -xqp; 0.0]; Cvec = Cvec/norm(Cvec);
        Rvec = [xqp;  yqp; 0.0]; Rvec = Rvec/norm(Rvec);
        
        fiber = Microstructure(thetaEndo, thetaMid, thetaEpi, Rendo, Repi, rqp, Rvec, Cvec, DebugFlag);
        
        [ThetaPol, Rpol, Zpol] = cart2pol(xqp, yqp, zqp);
        [xa, ya, za, Fan] = AnalyticalF(ThetaPol, Rpol, Zpol, Zbot, Ztop, afinal*alphaTime(s), betamax*alphaTime(s), OptBeta);
        [CellScalars(e, 1:NumInv), InvariantNames] = ComputeInvariants(Fan, fiber, Lvec, Cvec, Rvec);
        Eff_an(s,e) = CellScalars(e, 5);
        J_an(s,e) = sqrt(CellScalars(e, 3));
        FibersUpdated(e,:) = fiber(1:3)*transpose(Fan);
        FibersUpdated(e,:) = FibersUpdated(e,:)./norm(FibersUpdated(e,:));
        
        phiEl = xt(conn(e,:), :);
        
        numQP = 0;
        Fnum = LinHexa(Xel, phiEl, numQP);
        [CellScalars(e, NumInv+1:NumInv*2), InvariantNames] = ComputeInvariants(Fnum, fiber, Lvec, Cvec, Rvec);
        Eff_num(s,e) = CellScalars(e, NumInv+5);
        J_num(s,e) = sqrt(CellScalars(e, NumInv+3));
        
    end
    
    % Plot deformed heart
    ElType = 32;
    for (n = 1:NumInv)
        ScalarNames(n,:)        = [InvariantNames(n,:), '_An_'];
        ScalarNames(n+NumInv,:) = [InvariantNames(n,:), '_Num'];
    end
    PlotToVTK(xt, conn-1, ElType, [], [], CellScalars, ScalarNames, [OutputPathVTK,'PhantomAn'], s);
%     VectorsToVTK(nodesQP, [Eff_an(s,:)', Eff_num(s,:)'] , ['I4_an ';'I4_num'], FibersUpdated, 'Myofiber', [OutputPathVTK,'PhantomAnVec'])
end

figure(10);
hold on;
EffAnMedian = median(Eff_an,2);
plot(median(Eff_an,2) , '-or');
plot(median(Eff_num,2), '-ok');
legend('Analytical', 'Numerical');
xlabel('Time [ms]');
ylabel('E_{ff}');
hold off; 

figure(11);
hold on;
plot(median(J_an,2), '-or');
plot(median(J_num,2), '-ok');
legend('Analytical', 'Numerical');
xlabel('Time [ms]');
ylabel('J');
hold off;









