function result = func_U(s, omega, gamma, nls_case)
    a1 = nls_case.coefs(1);
    a3 = nls_case.coefs(2);
    if isequal(nls_case.powers,[3,5,7])
        result = omega*s-a1/2*s.^2+gamma/3*s.^3-a3/4*s.^4; % 357 case
    elseif isequal(nls_case.powers,[2,3,4])
        result = omega*s-2*a1/3*s.^(3/2)+gamma/2*s.^2-2*a3/5*s.^(5/2); % 234 case
    elseif isequal(nls_case.powers, [5,7,9])
        result = omega*s-a1/3*s.^3+gamma/4*s.^4-a3/5*s.^5; % 579 case
    end
end