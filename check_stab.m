clf;
clear;

nls_case.powers = [5,7,9];
nls_case.coefs = [1,-1]; % coefs of f; [a_1,a_3] in the notation of 234nls
log_scale = 1; % use log scale for omega; otherwise use linear scale
ignore_inf = 0; % replace infinite J values with NaN
tol = 0.00000001;

gamma_vals = linspace(-10,10,50);
if(log_scale)
    omega_vals = logspace(-3,1,50);
else
    omega_vals = linspace(0,5,50);
end
J_vals = zeros(length(gamma_vals),length(omega_vals));

num_issues = zeros(1,3); % number of (omega, gamma) such that [no a value found, Up is positive at a, U is negavite somewhere in (0,a)]
num_exist = 0; % number of (omega,gamma) such that a solution exists

idx_bdry = 1;
idx_stab = 1;
idx_unstab = 1;
idx_inf = 1;
idx_weird = 1;

stab_om=[];
stab_gam=[];
unstab_om=[];
unstab_gam=[];
boundary_om=[];
boundary_gam=[];


for i = 1: length(omega_vals)
    for j = 1: length(gamma_vals)
        J_vals(j,i) = NaN;
        om = omega_vals(i);
        gam = gamma_vals(j);
        alp = func_a(om,gam,nls_case);
        issues = DNE_omgam(om,gam,alp,nls_case,tol);
        if (~issues)
            Jval = func_J(om,gam,alp,nls_case);
            if abs(Jval) <= 10^(-6)
                boundary_om(idx_bdry) = om;
                boundary_gam(idx_bdry) = gam;
                idx_bdry = idx_bdry+1;
            elseif isinf(Jval)
                inf_om(idx_inf) = om;
                inf_gam(idx_inf) = gam;
                %disp(['Infinite J value for (om,gam)=', num2str(om), ' ', num2str(gam)])
                idx_inf = idx_inf+1;
                if ignore_inf
                    Jval = NaN;
                end
            elseif Jval>0
                stab_om(idx_stab) = om;
                stab_gam(idx_stab) = gam;
                idx_stab = idx_stab+1;
            elseif Jval < 0
                unstab_om(idx_unstab) = om;
                unstab_gam(idx_unstab) = gam;
                idx_unstab = idx_unstab+1;
            else
                weird_om(idx_weird) = om;
                weird_gam(idx_weird) = gam;
                weird_J(idx_weird) = Jval;
                idx_weird = idx_weird+1;
            end
            J_vals(j,i) = Jval;
            num_exist=num_exist+1;
        else
            num_issues(issues) = num_issues(issues) + 1;
        end
    end
end

if nls_case.coefs(1) == 1
    first_letter = 'F';
else
    first_letter = 'D';
end
if nls_case.coefs(2) == 1
    second_letter = 'F';
else
    second_letter = 'D';
end

hold on;

if log_scale
    plot(log10(stab_om),stab_gam,'o');
    plot(log10(unstab_om),unstab_gam,'*');
    %plot(log10(boundary_om),boundary_gam,'c');
    %plot(log10(inf_om),inf_gam,'x', 'Color','g');
    contour(log10(omega_vals),gamma_vals,J_vals,[-100,-8,-5,-2,-1,0,1,8,100],'showtext','on')
    xlabel('log10(omega)')
else
    plot(stab_om,stab_gam,'o');
    plot(unstab_om,unstab_gam,'*');
    plot(boundary_om,boundary_gam,'c');
    %plot(inf_om,inf_gam,'x', 'Color','g');
    contour(omega_vals,gamma_vals,J_vals,[-100,-9,-8,-1,-.1,-.01,0.01,.1,1,8,9,100],'showtext','on')
    xlabel('omega')
end
ylabel('gamma')
title(strcat(first_letter,second_letter, ' case: Curves and Regions'))

disp(['no a found: ', num2str(num_issues(1))])
disp(['Upa positive: ', num2str(num_issues(2))])
disp(['min U negative: ', num2str(num_issues(3))])


  

%{
for i=1: length(inf_om)
    inf_a(i) = func_a(inf_om(i),inf_gam(i),a1,a3);
end
display('Infinite J: (omega,gamma,a)')
display([inf_om;inf_gam;inf_a])
%}

%{
function result = J_integrand(s,omega,gamma,a,a1,a3)
    Uas = func_U(a*s,omega,gamma,a1,a3);
    Upa = func_Up(a,omega,gamma,a1,a3);
    Upas = func_Up(a*s,omega,gamma,a1,a3);
    result = (3+a*s.*(Upa-Upas)/Uas).*realsqrt(s)/realsqrt(Uas);
end
%}
%{
function result = func_U(s, omega, gamma, a1, a3)
    result = 2*func_GG_sqrt(s,omega,gamma,a1,a3);
end

function result = func_Up(s,omega,gamma,a1,a3)
    result = omega-a1*s+gamma*s.^2-a3*s.^3;
end
%}

