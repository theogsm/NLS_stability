function display_Jgam(om,gamBounds, a1,a3)
    J = @(g) func_J(om,g,func_a(om,g,a1,a3),a1,a3);
    interval = linspace(gamBounds(1),gamBounds(2),100);
    plot(interval,arrayfun(J,interval));
end