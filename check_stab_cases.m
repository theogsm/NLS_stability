%function [iter, stab_ome, stab_gam,instab_ome,instab_gam,infnorm_ome,infnorm_gam,bound_weird_ome,bound_weird_gam,boundary_ome,boundary_gam,upa_region,ua_region] = check_stab_cases(omega_exist,gamma_exist,b,d)
%% use this code for case 3, always exist, but numerically, we see that there exists stability boundary curve.
gam1 = 4/(sqrt(5));
ome1 = 2*(sqrt(5))/27;
gamma_exist = linspace(-5,5,50);
omega_exist = linspace (0,5,50);

gamma_002 = gamma_exist;
omega_002 = omega_exist;

format long g;
acc = 1;
looping = 1;
acc_b= 1;
acc_stab = 1;
acc_instab = 1;
acc_inf = 1;
no_A_acc = 1;
norm_acc_weird =1;
acc_weird = 1; 
MaxQuadGK = 1200;

b = -1;
d = 1;

% % plot the non-exist curve:
% alp_b = sqrt(5/3);
% %alp = linspace(alp_b,10, 100);
% alp = linspace(0,alp_b, 100);


% gamma = -(2/3).*(alp.^-1)-(6/5).*alp;
% omega = -(1/3).*alp + (alp.^(3))/5;
% 
% alp_crit = alp_b;
% gamma_crit = -(2/3).*(alp_crit.^-1)-(6/5).*alp_crit;
% omega_crit = -(1/3).*alp_crit + (alp_crit.^(3))/5;
%plot(omega, gamma, 'black');
%plot(omega_crit,gamma_crit, 'o');
% xlabel('omega');
% ylabel('gamma');
% title('Nonexist curve, 0<a<5/9');

% using c = -gamma, we plot omega vs c for non-exist
%c = (2/3).*(alp.^-1)+(6/5).*alp;
%plot(omega,c)

tol = 0.00000001;
acci = 1;
accj = 1;

% initializing matrices:
norm_N = length(omega_exist);
a_exist = -1994.668*ones(norm_N);
norm_exist = 1994.668*ones(norm_N);

stab_ome = [];
stab_gam = [];
instab_ome = [];
instab_gam = [];
norm_exist = [];

boundary_ome = [];
boundary_gam = [];


for i = 1: length(omega_exist)
    for j = 1: length (gamma_exist)
        disp('start looping to find existing omega and gamma...');
        try func_a_cases(omega_exist(i),gamma_exist(j),1,b,d)
            alp = func_a_cases(omega_exist(i),gamma_exist(j),1,b,d);
        upa_region_temp = omega_exist(i) + func_f_cases(alp,gamma_exist(j),b,d);
        disp('a root found');
        ua_region(i,j) = func_u_cases(alp,omega_exist(i),gamma_exist(j),b,d);
        upa_region_matrix(j,i) = upa_region_temp;
        epsii = (alp-0)/100;
        s = linspace(0+epsii,alp-epsii,50);
        val = func_u_cases(s,omega_exist(i),gamma_exist(j),b,d);
        min_us_temp = min(val);
        min_us_matrix(j,i) = min_us_temp;

        if upa_region_temp <= tol
%             s = linspace(0+0.05,alp-0.05,50); % epsilon to nonexist bound
            disp('upa<0 checked');
%             val = func_u_cases(s,omega_exist(i),gamma_exist(j),b,d);
%             min_us_temp = min(val);
            if (min_us_temp >= tol)% u(s) >0, and intro acc for existence ome,gam
            exist_curve_ome(acc) = omega_exist(i);
            exist_curve_gam(acc) = gamma_exist(j);
            
            if abs(alp) <= tol
               integrand = @(s) 3+s^(-0.5).*(1-s)./(2/3) % limit of the integrand when a <<1
            
            else
                integrand = @(s) dQdw_integrand_cases(s,omega_exist(i),gamma_exist(j),alp,b,d);  
            end% integrand, to be taken for s in [0,1]
            norm_exist_temp = -0.5*alp^1.5/upa_region_temp*quadgk(integrand,0,1,'MaxIntervalCount',MaxQuadGK);
             
            
            if abs(norm_exist_temp) <= 10^(-6)  %% for dqdw close  zero
                boundary_ome(acc_b) = omega_exist(i);
                boundary_gam(acc_b) = gamma_exist(j);
                acc_b = acc_b +1;
                
            elseif isinf(norm_exist_temp)
                infnorm_ome(acc_inf) = omega_exist(i);
                infnorm_gam(acc_inf) = gamma_exist(j);
                acc_inf = 1 + acc_inf;
                
            elseif norm_exist_temp >0
                    stab_ome(acc_stab) = omega_exist(i);
                    stab_gam(acc_stab) = gamma_exist(j);
                    acc_stab = acc_stab +1;
            elseif norm_exist_temp <0
                    instab_ome(acc_instab) = omega_exist(i);
                    instab_gam(acc_instab) = gamma_exist(j);
                    acc_instab = acc_instab +1;
                    
            else 
                bound_weird_ome(acc_weird) = omega_exist(i);
                bound_weird_gam(acc_weird) = gamma_exist(j);
                norm_weird (norm_acc_weird) = norm_exist(acci,accj);
                norm_acc_weird = 1 + norm_acc_weird;
                acc_weird = 1+ acc_weird;
            end
            
           acc = acc+1;
           disp('u(s)>0 checked');
            else
                problem_us(problem_us_acc,1) = (omega_exist(i));
                problem_us(problem_us_acc,2) = (gamma_exist(h));
                problem_us_acc = 1+ problem_us_acc;
                acc = acc +1;
            end
        
            
        end
        norm_exist(acci,accj) = norm_exist_temp;
        a_exist(acci,accj) = alp;
        acci = 1+ acci;
    
        
        catch
            disp('no a found');
            no_A_acc = 1+ acc;
        end 
        
    looping = looping + 1
    end
    acci = 1;
    accj = 1 + accj;
end

iter = [looping,acc_b,acc_stab,acc_instab,acc_inf,norm_acc_weird,acc_weird]
% 
% store fine:
% stab_ome_fine = stab_ome;
% stab_gam_fine = stab_gam;
% instab_ome_fine = instab_ome;
% instab_gam_fine = instab_gam;
% norm_fine = norm_exist;
% bound_ome_fine = boundary_ome;
% bound_gam_fine = boundary_gam;



% store coarse:
stab_ome_002 = stab_ome;
stab_gam_002 = stab_gam;
instab_ome_002 = instab_ome;
instab_gam_002 = instab_gam;
norm_002 = norm_exist;

bound_ome_002 = boundary_ome;
bound_gam_002 = boundary_gam;

figure
hold on
plot(stab_ome,stab_gam,'o')
plot(instab_ome, instab_gam, '*')
plot(boundary_ome,boundary_gam,'c');

%plot(infnorm_ome, infnorm_gam,'s')
%title('D*F case Stab regions;exists cond 1.5-1.6')
% xlabel('omega')
% ylabel('gamma')
%legend('stable','instable','numerical boundary','dQdw-level curves');

% Plot the pole for our case 
%plot([2*sqrt(5)/27],[4/sqrt(5)],'kx','markersize',10)
%plot(omega, gamma, 'r');
%end

% level curves of dNorm:
contour(omega_exist,gamma_exist,norm_exist,[-100,-9,-8,-0.4,0,0.05,0.1,0.4,8,9,100],'showText','on')
title('D*F case Stab regions;exists cond 1.5-1.6 and dQdw contour plot')
xlabel('omega')
ylabel('gamma')
text(20,-2,'Sigma_S','FontSize',14,'Color','blue')
text(20,4,'Sigma_U','FontSize',14,'Color','red')
legend('stable','instable','numerical boundary','dQdw-level curves');


function t = func_u_cases(s,omega,gamma,b,d)
 t = omega.*s-(2/3)*b.*s.^(1.5)+(gamma/2).*s.^2-(2/5)*d.*s.^2.5;
 t = real(t);
end

%
% Function : dQdw
%
% compute dQdw for a given value of omega,gamma; expected as a 2-vector 'input'
%    
function [func,gradient] = dQdw_func_cases(input,b,d)
    global DEBUG; MaxQuadGK = 1200;
    omega = input(1);
    gamma = input(2);
    if(DEBUG) disp('BEGIN dQdw:'); disp('Calling computeA(...)'); end
    a = func_a(omega,gamma);
    Upa = omega + func_f_cases(a,gamma,b,d);
    
    if(DEBUG) disp('Calling quadgk(...)'); end
    integrand = @(s) dQdw_integrand_cases(s,omega,gamma,a,b,d);      % integrand, to be taken for s in [0,1]
    func = -0.5*a^1.5/Upa*quadgk(integrand,0,1,'MaxIntervalCount',MaxQuadGK);      % see 'Explicit Formulas', page 1

    if nargout > 1
        gradient = zeros(1,2);
        if(DEBUG) disp('Calling quadgk(...) for dQdww'); end
        integrand = @(s) dQdww_integrand_cases(s,omega,gamma,a,b,d);      % integrand, to be taken for s in [0,1]       
        gradient(1) = -0.5*a^0.5/Upa*func ...
                       +0.5*a^1.5/Upa^2*quadgk(integrand,0,1,'MaxIntervalCount',MaxQuadGK);  
                       % see 'Explicit Formulas', page 2          
        if(DEBUG) disp('Calling quadgk(...) for dQdgw'); end       
        integrand = @(s) dQdgw_integrand_cases(s,omega,gamma,a,b,d);      % integrand, to be taken for s in [0,1]       
        gradient(2) = -0.5*func_Gg(a)/a^0.5/Upa*func ...
                       +0.5*a^1.5/Upa^2*quadgk(integrand,0,1,'MaxIntervalCount',MaxQuadGK);  
                       % see 'Explicit Formulas', page 4          
    end
    
    if(DEBUG) disp('END dQdw'); end
end % of function dQdw

function t = func_f_cases(s,gamma,b,d) %(f1)
    t = -b*s.^0.5 + gamma*s - d*s.^1.5;  %2.1 equ
end % of user-supplid function func_f(s)

function result = dQdw_integrand_cases(s,omega,gamma,a,b,d)
    as = a*s;
    Uas = arrayfun(@(x) max(omega*x+func_g_cases(x,gamma,b,d),0), as);
    UpaUpas = func_f_cases(a,gamma,b,d)- arrayfun(@(x) func_f_cases(x,gamma,b,d),as);
    result = ( 3 + as.*UpaUpas./Uas ) ...
             .* realsqrt(s) ./ realsqrt(Uas);   
end % of function dQdw_integrand

function t = func_g_cases(s,gamma,b,d)
    t = -(2/3)*b.*s.^1.5 + gamma/2.*s.^2 - (2/5)*d.*s.^2.5;
end % of user-supplid function func_g(s)

function result = dQdww_integrand_cases(s,omega,gamma,a,b,d)
    as = a*s;
    Uas = arrayfun(@(x) max(omega*x + func_g_cases(x,gamma,b,d),0), as); 
    Upa = omega + func_f_cases(a,gamma,b,d);
    Upas = arrayfun(@(x) omega + func_f_cases(x,gamma,b,d), as);
    UpaUpas = Upa - Upas;
    Uppa = func_df(a,gamma); % derivative of f(s,gamma) by s
    Uppas = arrayfun(@(x) func_df_cases(x,gamma,b,d),as);
     
    result = (   3*(1 - a*Uppa/Upa) + ...
                 as.*UpaUpas./Uas.*( 2.5 + 1.5*as.*UpaUpas./Uas ) + ...
                 -as.*( (as.*Uppas-Upa)+(Upa-a*Uppa)*Upas/Upa )./Uas ...
                 ) ...
             .* realsqrt(s) ./ realsqrt(Uas);   
end % of function dQdww

function t = func_df_cases(s,gamma,b,d) % derivative of f(s,gamma) by s
    tau=s^0.5;
    t = -0.5*b/tau + gamma - d*1.5*tau;
end % of user-supplid function func_df(s)    

function result = dQdgw_integrand_cases(s,omega,gamma,a,b,d)
    as = a*s;
    Uas = arrayfun(@(x) max(omega*x + func_g_cases(x,gamma,b,d),0), as); 
    Upa = omega + func_f_cases(a,gamma,b,d);
    Upas = arrayfun(@(x) omega + func_f_cases(x,gamma,b,d), as);
    UpaUpas = Upa - Upas;
    Uppa = func_df_cases(a,gamma,b,d); % derivative of f(s,gamma) by s
    Uppas = arrayfun(@(x) func_df_cases(x,gamma,b,d),as);

    Gga = func_Gg(a,gamma);  % derivative of g(s,gamma) by gamma, nothing involved with b, d. No changes 
    Ggas = arrayfun(@(x) func_Gg(x,gamma),as);
    Gfa = func_Gf(a,gamma); % derivative of f(s,gamma) by gamma = derivative of g(s,gamma) by s and gamma
    Gfas = arrayfun(@(x) func_Gf(x,gamma),as);
    
    result = (   3*Gga*(1-a*Uppa/Upa) + 3*Gfa*a ...
                 + 1.5*a*(Ggas*Upa - Gga*s.*Upas)./Uas ...
                 + as.*UpaUpas./Uas*Gga ...
                 + 1.5*a*(Ggas*Upa - Gga*s.*Upas)./Uas.*as.*UpaUpas./Uas ...
                 + s.*UpaUpas./Uas*Gga ...
                 + a*Gga*s.*s.*(Uppa/Upa*Upas - Uppas)./Uas ...
                 + as.*(Upa*Gfas-Upas*Gfa)./Uas ...
                 ) ...
             .* realsqrt(s) ./ realsqrt(Uas);      
end % of function dQdgw
    

function a = func_a_cases(omega,gamma,varargin,b,d)

    % Seek 0 = U(a) = U(a)/a = omega + g(a)/a
    % Note that g(a) is a cubic polynomial in a^0.5
    possibleRoots = cubicroots( +(2.0/5.0)*d, -gamma/2.0, +(2.0/3.0)*b, -omega);

    % The roots returned are ALWAYS real-valued; take the smallest positive one, and square it
    possibleRoots = sort(possibleRoots);
    idx = min(find(possibleRoots > 0));
    if(isempty(idx)) error('No root found; cannot compute a'); end
    a = ( possibleRoots(idx) )^2;
    
% $$$     persistent init;
% $$$     if(~isempty(varargin)) init = varargin{1};
% $$$     elseif(isempty(init)) error('Need to initialize before use');
% $$$     end % else, do nothing, use saved value
% $$$     a = fsolve(@(x) omega*x + g(x,gamma),init); 
% $$$     if(~isreal(a)) error('Complex result indicates poor initial point'); end
end 