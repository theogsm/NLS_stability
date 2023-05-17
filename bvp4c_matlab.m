close all;
clear;
clc;
% Setting omega, p
om=1;
p=3;

% Set u0 to be the first zero of G
u0=((p+1)*om/2).^(1/(p-1));

xmesh = linspace(0,10,5);

solinit = bvpinit(xmesh, @guess);

% use bvp4c for the solution (u,u') passing through (u0,0)
sol = bvp5c(@(x,y)(bvpfcn(x,y,om,p)), @(ya,yb)(bcfcn(ya,yb,u0)), solinit);

plot(sol.x, sol.y(1,:))
hold on;
exact_y=exact_sol(sol.x,om,p);
plot(sol.x, exact_y)
hold off;

%plot(sol.x, sol.y(1,:)-exact_y);

% functions encoding the equation an boundary conditions for passing to bvp4c
function dydx = bvpfcn(x,y,om,p)
    dydx = zeros(2,1);
    dydx = [y(2);om*y(1)-y(1).^p];
end

function res = bcfcn(ya,yb, u0)
    res = [ya(1)-u0;ya(2)];
end

% solution initial guess
function g = guess(x)
    g = [sin(x);cos(x)];
end

% exact solution given in section 6 of 234nls
function u = exact_sol(x,om,p)
    u=om.^(1/(p-1))*Q(om.^(1/2)*x,p);
end

function q = Q(x,p)
    q=((p+1)/2*(sech((p-1)*x/2)).^2).^(1/(p-1));
end
