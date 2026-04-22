clearvars
close all
clc 
% Inputs
n     = 31;        % No. of grid points       
nt    = 1400;       % No. of time steps
gamma = 1.4;        % Heat capacity ratio
C     = 0.5;        % Courant Number 
tol   = 1e-8;       % Error tolerance at convergence

% Solving 
[rho1,v1,T1,p1,M1,mflow_nc,kn1] = non_conservation_form(n,nt,C,gamma,tol);
[rho2,v2,T2,p2,M2,mflow_c,kn2] = conservation_form(n,nt,C,gamma,tol);

x = linspace(0,3,n);
str1 = ["(Non-Conservation form)" "(Conservation form)"];
 
for j = 1:2
 % Distributing variables   
 if  j==1 , [rho,v,T,p,M,mflow,kn] = deal(rho1,v1,T1,p1,M1,mflow_nc,kn1);
 else ,    [rho,v,T,p,M,mflow,kn] = deal(rho2,v2,T2,p2,M2,mflow_c,kn2);end
     
% BLOCK 1: Plotting Steady state distribution of variables inside nozzle
figure('position',[300 100 1000 400])
plot(x,v(1,:,kn),'b',x,T(1,:,kn),'r',x,rho(1,:,kn),'g',x,p(1,:,kn),'m')
legend({'Velocity','Density','Temperature','Pressure'},'location','northwest')
xlabel('Nozzle length ( x / L ) rightarrow')
ylabel('Non-dimensional variables rightarrow')
title('Steady state distribution of primitive variables')
subtitle(str1(j))


% BLOCK 2: Plotting timewise variation of massflow rate along nozzle length
figure('position',[300 100 1000 500])
clr = ["-.b" "c" "y" "g" "m" "xr" "k"];
i   = 1;
for k = 1:nt
    if ismember(k,[1 50 100 150 200 700 kn])
    plot(x,mflow(:,:,k),clr(i),'linewidth',0.9)
    hold on
    i = i+1;
    end
end
str2 = ["1Deltat","50Deltat","100Deltat",...
    "150Deltat","200Deltat","700Deltat",kn + "Deltat"];
legend(str2,'location','north','NumColumns',2)
xlabel('Nozzle length ( x / L ) rightarrow')
ylabel('Mass flow rate  ( rho / rho* ) rightarrow')
title('Transient variation of Mass flow rate (Non-dimensional)')
subtitle(str1(j));
hold off

end

% BLOCK 3: Plotting timewise variaton of variables at throat
% getting throat index
if (mod(n,2)==0),throat = n/2; % n is even
else ,throat = (n+1)/2;        % n is odd         
end
figure('position',[300 100 1000 500])
t = tiledlayout(4,1);

nexttile
plot(1:kn1, squeeze(v1(1,throat,:)) ,'g','linewidth',0.8);
hold on
plot(1:kn2, squeeze(v2(1,throat,:)) ,'k','linewidth',0.8);
title('Velocity')
legend('Non-conservation form','Conservation form')

nexttile
plot(1:kn1, squeeze(p1(1,throat,:)) ,'r','linewidth',0.8);
hold on
plot(1:kn2, squeeze(p2(1,throat,:)) ,'k','linewidth',0.8);
title('Pressure')
legend('Non-conservation form','Conservation form')

nexttile
plot(1:kn1, squeeze(T1(1,throat,:)) ,'b','linewidth',0.8);
hold on
plot(1:kn2, squeeze(T2(1,throat,:)) ,'k','linewidth',0.8);
title('Temperature')
legend('Non-conservation form','Conservation form')

nexttile
plot(1:kn1, squeeze(rho1(1,throat,:)) ,'m','linewidth',0.8);
hold on
plot(1:kn2, squeeze(rho2(1,throat,:)) ,'k','linewidth',0.8);
title('Density')
legend('Non-conservation form','Conservation form')

xlabel(t,'Time steps rightarrow')
ylabel(t,'Non-dimensional variable rightarrow')
title(t,'Transient variation of Non-dimensional variables')
subtitle(t,'(At nozzle throat)')

%BLOCK 4: Comparison of Normalized mass flow rate distributions of both forms
figure('position',[300 100 1000 500]) 
plot(x,mflow_nc(1,:,kn1),'m','linewidth',1)
hold on
plot(x,mflow_c(1,:,kn2),'b','linewidth',1)
plot(x,0.579*ones(1,n),'g','linewidth',1)
legend({'Non-Consevation Form','Consevation Form','Analytical Solution'},'location','north')
title('Comparison of Normalized mass flow rate distributions of both forms')
subtitle('(Steady state)')
xlabel('Nozzle length ( x / L ) rightarrow')
ylabel('Mass flow rate  ( rho / rho* ) rightarrow')

row_tab1 = ["Non-conservation form";"Conservation form"];
Steady_state_mean_mfr  = [mean(mflow_nc(1,:,kn1));mean(mflow_c(1,:,kn2))];
Analytical_mfr         = [0.579;0.579];
Error_percent          = 100*(Steady_state_mean_mfr - Analytical_mfr);
Iterations_to_converge = [kn1;kn2];
tab1 = table(Steady_state_mean_mfr,Analytical_mfr,Error_percent,...
                Iterations_to_converge,'rownames',row_tab1);
disp(tab1)

% BLOCK 5: Grid Independence test
gp = [31,62,93];
for j = 1:3
% [rho_t,v_t,T_t,p_t,M_t,mflow_t,kn]
[rho1,v1,T1,p1,M1,mfr1,kn1]= non_conservation_form(gp(j),nt,C,gamma,tol);
[rho2,v2,T2,p2,M2,mfr2,kn2]  = conservation_form(gp(j),nt,C,gamma,tol);
% getting throat index
if (mod(gp(j),2)==0),throat = gp(j)/2; % n is even
else ,throat = (gp(j)+1)/2; end        % n is odd         

Density_NC(j,1)     = rho1(1,throat,kn1);
Temperature_NC(j,1) = T1(1,throat,kn1);
Pressure_NC(j,1)    = p1(1,throat,kn1);
Mach_no_NC(j,1)     = M1(1,throat,kn1);

Density_C(j,1)     = rho2(1,throat,kn2);
Temperature_C(j,1) = T2(1,throat,kn2);
Pressure_C(j,1)    = p2(1,throat,kn2);
Mach_no_C(j,1)     = M2(1,throat,kn2);
end

fprintf('n                              [GRID INDEPENDENCE TEST]')
fprintf('nFor Non-conservation form:n')
Grid_Points = num2str(gp');
tab2 = table(Grid_Points,Density_NC,Temperature_NC,Pressure_NC,Mach_no_NC);
disp(tab2)
fprintf('For Conservation form:n')
tab3 = table(Grid_Points,Density_C,Temperature_C,Pressure_C,Mach_no_C);
disp(tab3)
% END OF PROGRAM