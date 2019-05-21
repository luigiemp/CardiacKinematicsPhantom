function G = EulerianMapNoOptToolbox(X, xgiven, Zbot, Ztop, a_s, beta_s, OptBeta)

[Phi, R, Z] = cart2pol(X(1), X(2), X(3));

[xa, ya, za, Fan] = AnalyticalF(Phi, R, Z, Zbot, Ztop, a_s, beta_s, OptBeta);

G = norm(xgiven - [xa, ya, za]);

end

