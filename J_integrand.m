function result = J_integrand(s,omega,gamma,a,nls_case) % integrand from (2.12) in 234nls
    Uas = arrayfun(@(s)max(func_U(s,omega,gamma,nls_case),0), a*s); % U(as)
    Upa = func_Up(a,omega,gamma,nls_case); % U'(a)
    Upas = arrayfun(@(s) func_Up(s,omega,gamma,nls_case),a*s); % U'(as)
    result = (3+a*s.*(Upa-Upas)./Uas).*realsqrt(s)./realsqrt(Uas);
end