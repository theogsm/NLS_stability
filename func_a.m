function a = func_a(omega, gamma, a1, a3)
    powers = '357';
    if strcmp(powers,'357')
        roots = cubicroots(-a3/4, gamma/3, -a1/2, omega); % 357 case
    elseif strcmp(powers,'234')
        roots = cubicroots( 2/5*a3, -gamma/2, 2/3*a1, -omega); % 234 case
    end

    smallest = min(roots(roots>0));

    if(isempty(smallest))
        a=-1;
    else
        if strcmp(powers,'357')
            a = smallest; % 357 case
        elseif strcmp(powers,'234')
            a = smallest^2; % 234 case
        end
        
    end
end
