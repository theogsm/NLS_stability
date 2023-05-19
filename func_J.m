function result = func_J(omega,gamma,alp,a1,a3) % fromula for d'' from (2.12) in 234nls
    MaxIntervals = 10000;
    integrand = @(s) J_integrand(s,omega,gamma,alp,a1,a3);
    Upa = func_Up(alp,omega,gamma,a1,a3);
    result = -0.5*alp^1.5/Upa*quadgk(integrand,0,1,"MaxIntervalCount",MaxIntervals);
end