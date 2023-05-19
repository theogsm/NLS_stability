function result = func_U(s, omega, gamma, a1, a3)
    powers = '357';
    if strcmp(powers,'357')
        result = omega*s-a1/2*s.^2+gamma/3*s.^3-a3/4*s.^4; % 357 case
    elseif strcmp(powers,'234')
        result = omega*s-2*a1/3*s.^(3/2)+gamma/2*s.^2-2*a3/5*s.^(5/2); % 234 case
    end
end