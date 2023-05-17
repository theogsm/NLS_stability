function result = special_func_J(om, gam, alp, a1, a3)
    % only for 357 case
    MaxIntervals =50000;
    integrand = @(u) arrayfun(@(s) special_J_integrand(s, om,gam,alp,a1,a3),u);
    Upa = func_Up(alp,om,gam,a1,a3);
    result = -0.5*alp^1.5/Upa*quadgk(@(u) arrayfun(integrand,u),1,Inf,"MaxIntervalCount",MaxIntervals);
end