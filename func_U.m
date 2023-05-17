function result = func_U(s, omega, gamma, a1, a3)
    % 357 case:
    result = omega*s-a1/2*s.^2+gamma/3*s.^3-a3/4*s.^4;
    
    % 234 case:
    %result = omega*s-2*a1/3*s.^(3/2)+gamma/2*s.^2-2*a3/5*s.^(5/2);
end