function result = func_Up(s,omega,gamma,a1,a3)
    powers = '357';
    if strcmp(powers,'357')
        result = omega - a1*s + gamma*s.^2 - a3*s.^3; % 357 case
    elseif strcmp(powers,'234')
        result = omega - a1*s.^.5 + gamma*s - a3*s.^1.5; % 234 case
    end
end