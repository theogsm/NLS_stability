function dne = DNE_omgam(om,gam,alp,a1,a3,tol)
    if alp >-.5 % func_a returns -1 if no root is found
        %disp('root of U found')
        Upa = func_Up(alp,om,gam,a1,a3);
        Upz = func_Up(0,om,gam,a1,a3);
        
        %eps = min(.5*alp,max(100*tol/Upa,100*tol/Upz));
        eps = alp/100;
        interval = linspace(eps,alp-eps,50);
        Uinterval = func_U(interval,om,gam,a1,a3);
        minU = min(Uinterval);
        if Upa < 0
            %disp('Upa < 0 checked')
            if minU > 0
                %disp('U>0 on (0,alp) checked')
                dne = 0; % solution exists
            else
                dne = 3; % solution DNE; U>0 fails
                disp(minU)
            end
        else
            dne = 2; % solution DNE; Upa<0 fails
        end
    else
        dne = 1; % solution DNE; no a value found
    end
end