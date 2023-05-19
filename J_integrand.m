function result = J_integrand(s,omega,gamma,a,a1,a3) % integrand from (2.12) in 234nls
    Uas = arrayfun(@(s)max(func_U(s,omega,gamma,a1,a3),0), a*s); % U(as)
    Upa = func_Up(a,omega,gamma,a1,a3); % U'(a)
    Upas = arrayfun(@(s) func_Up(s,omega,gamma,a1,a3),a*s); % U'(as)
    result = (3+a*s.*(Upa-Upas)./Uas).*realsqrt(s)./realsqrt(Uas);
end