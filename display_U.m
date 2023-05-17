function display_U(om, gam, a1,a3)
    clf;
    alp = func_a(om,gam,a1,a3);
    U = @(s) func_U(s,om,gam,a1,a3);
    interval = linspace(-alp/2,1.5*alp,100);
    hold on
    grid on
    plot(interval, U(interval));
    plot(alp,0,'o')
    
end