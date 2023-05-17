function result = special_J_integrand(u,om,gam,alp, a1,a3)
    s=@(v)(1-1/v);
    denom = max(om*alp-a1/2*alp^2*s(u)+gam/3*alp^3*s(u)^2-a3/4*alp^4*s(u)^3,0);
    % if denom<0
    %     disp('negative radicand in J integrand')
    %     [om,gam,alp,denom]
    % end
    result = (3+alp*(-a1*alp+gam*alp^2-a3*alp^3+a1*alp*s(u)-gam*alp^2*s(u)^2+a3*alp^3*s(u)^3)/(om*alp-a1/2*alp^2*s(u)+gam/3*alp^3*s(u)^2-a3/4*alp^4*s(u)^3))*(realsqrt(denom))^(-1)*u^(-2);
end