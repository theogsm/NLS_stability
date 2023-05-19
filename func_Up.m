function result = func_Up(s,omega,gamma,nls_case)
    a1 = nls_case.coefs(1);
    a3 = nls_case.coefs(2);
    if isequal(nls_case.powers,[3,5,7])
        result = omega - a1*s + gamma*s.^2 - a3*s.^3; % 357 case
    elseif isequal(nls_case.powers,[2,3,4])
        result = omega - a1*s.^.5 + gamma*s - a3*s.^1.5; % 234 case
    elseif isequal(nls_case.powers,[5,7,9])
        result = omega - a1*s^2 + gamma*s.^3 - a3*s.^4; % 579 case
    end
end