function a = func_a(omega, gamma, nls_case)
    a1 = nls_case.coefs(1);
    a3 = nls_case.coefs(2);
    if isequal(nls_case.powers, [3,5,7])
        Uroots = cubicroots(-a3/4, gamma/3, -a1/2, omega); % 357 case
    elseif isequal(nls_case.powers, [2,3,4])
        Uroots = cubicroots( 2/5*a3, -gamma/2, 2/3*a1, -omega); % 234 case
    elseif isequal(nls_case.powers, [5,7,9])
        Uroots = roots([-a3/5;gamma/4;-a1/3;0;omega]);
        Uroots = Uroots(imag(Uroots) == 0);
    end

    smallest = min(Uroots(Uroots>0));

    if(isempty(smallest))
        a=-1;
    else
        if isequal(nls_case.powers, [3,5,7]) || isequal(nls_case.powers, [5,7,9])
            a = smallest; % 357 case
        elseif isequal(nls_case.powers, [2,3,4])
            a = smallest^2; % 234 case
        end
        
    end
end
