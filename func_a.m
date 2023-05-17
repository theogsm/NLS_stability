function a = func_a(omega, gamma, a1, a3)
% 357 case:
roots = cubicroots(-a3/4, gamma/3, -a1/2, omega);

% 234 case:
%roots = cubicroots( 2/5*a3, -gamma/2, 2/3*a1, -omega);

smallest = min(roots(roots>0));

if(isempty(smallest))
    a=-1;
else
    % 357 case:
    a = smallest;

    % 234 case:
    %a = smallest^2;
end
end
