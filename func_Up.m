function result = func_Up(s,omega,gamma,a1,a3)
    % 357 case:
    result = omega - a1*s + gamma*s.^2 - a3*s.^3;
    
    % 234 case:
    %result = omega - a1*s.^.5 + gamma*s - a3*s.^1.5;
end