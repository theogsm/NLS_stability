function result = J_integrand(s,omega,gamma,a,a1,a3)
    Uas = arrayfun(@(s)max(func_U(s,omega,gamma,a1,a3),0), a*s);
    Upa = func_Up(a,omega,gamma,a1,a3);
    Upas = arrayfun(@(s) func_Up(s,omega,gamma,a1,a3),a*s);
    result = (3+a*s.*(Upa-Upas)./Uas).*realsqrt(s)./realsqrt(Uas);
end