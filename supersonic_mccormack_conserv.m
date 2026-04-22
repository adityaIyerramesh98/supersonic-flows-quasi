function [rho_t,v_t,T_t,p_t,M_t,mflow_t,kn] = conservation_form(n,nt,C,g,tol)

x = linspace(0,3,n);              % Space Grid
dx  = x(2) - x(1);                % Space-step
A   = 1 + 2.2*(x - 1.5).^2;       % Area profile (Non-Dimensional)
% getting throat index
if (mod(n,2)==0),throat = n/2; % n is even
else ,throat = (n+1)/2;        % n is odd         
end

% Preallocating arrays for speed
[du1dt,du2dt,du3dt,du1dt_p,du2dt_p,du3dt_p] = deal(zeros(1,n-1));
[rho_t(:,:,1),v_t(:,:,1),T_t(:,:,1),...
      p_t(:,:,1),M_t(:,:,1),mflow_t(:,:,1)] = deal(zeros(1,n));
                                    [rho,T] = deal(zeros(1,n));
                                    
% Initial Conditions
for z = 1:n
    if x(z)<=0.5
        rho(z) = 1;               % Density (Non-Dimensional)
        T(z)   = 1;               % Temprature (Non-Dimensional)
    elseif (x(z)>0.5) && (x(z)<=1.5)
        rho(z) = 1-0.366*(x(z) - 0.5);
        T(z)   = 1 - 0.167*(x(z) - 0.5);
    elseif (x(z)>1.5) && (x(z)<=3.5)
        rho(z) = 0.634 - 0.3879*(x(z) - 1.5);
        T(z)   = 0.833 + 0.3507*(x(z) - 1.5);
    end
end              
v   = (0.59)./(rho.*A); % Velocity (Non-Dimensional)
u1  = rho.*A;
u2  = 0.59*ones(1,n);
u3  = u1.*((T/(g-1)) + (g/2).*(v.^2));

% Conservation form of Governing Equ.
   for k = 1:nt
        dt  = min((C*dx)./((T.^0.5)+v));  % Selecting min value of dt 
                                          % across all grid points
        % Updating old values
        u1_old  = u1;
        u2_old  = u2;
        u3_old  = u3;
        % Flux terms at kth time level
        f1 = u2;
        f2 = ((u2.^2)./u1) +((g-1)/g)*(u3 - (g/2)*(u2.^2)./u1);
        f3 = ((g*u2.*u3)./u1) - ((g*(g-1)/2)*(u2.^3)./(u1.^2));

        % Predictor Step
        for i = 2:n-1
            % Foreward differences (at ith grid point)
            % Source term   
            dAdx = (A(i+1)-A(i))/dx; 
            j2   = ((g-1)/g)*(u3(i) - 0.5*g*((u2(i)^2)/u1(i)) )*dAdx/A(i);

            % Continuity equation
            du1dt(i) = -(f1(i+1) - f1(i))/dx;
            % Momentum equation
            du2dt(i) = -(f2(i+1) - f2(i))/dx + j2;
            % Energy equation
            du3dt(i) = -(f3(i+1) - f3(i))/dx;

            % Solution update (predicted values)
            u1(i) = u1_old(i) + du1dt(i)*dt;
            u2(i) = u2_old(i) + du2dt(i)*dt;
            u3(i) = u3_old(i) + du3dt(i)*dt;

            % Predicted values of fluxes
            f1(i) = u2(i);
            f2(i) = ((u2(i)^2)/u1(i)) +((g-1)/g)*(u3(i) - (g/2)*(u2(i)^2)/u1(i));
            f3(i) = ((g*u2(i)*u3(i))/u1(i)) - ((g*(g-1)/2)*(u2(i)^3)/(u1(i)^2));
        end

        % Corrector Step
        for j = 2:n-1
            % Backward differences (at jth grid point)       
            % Source term        
            dAdx = (A(j)-A(j-1))/dx;
            j2   = ((g-1)/g)*(u3(j) - 0.5*g*((u2(j)^2)/u1(j)))*dAdx/A(j);

            % Continuity equation
            du1dt_p(j) = -(f1(j) - f1(j-1))/dx;
            % Momentum equation
            du2dt_p(j) = -(f2(j) - f2(j-1))/dx + j2;
            % Energy equation
            du3dt_p(j) = -(f3(j) - f3(j-1))/dx;
        end

        % Average Time derivatives
        du1dt_av   = 0.5*(du1dt + du1dt_p); % drhodt_av
        du2dt_av   = 0.5*(du2dt + du2dt_p);
        du3dt_av   = 0.5*(du3dt + du3dt_p);

        % Final Values at (k+1)th time step 
        for m = 2:n-1
            u1(m) = u1_old(m) + du1dt_av(m)*dt;
            u2(m) = u2_old(m) + du2dt_av(m)*dt;
            u3(m) = u3_old(m) + du3dt_av(m)*dt;
        end

        % Boundary Conditions
        % Inlet boundary condition
        u1(1)  = A(1);
        u2(1)  = 2*u2(2) - u2(3);
        u3(1)  = u1(1)*( (T(1)/(g-1)) + (g/2)*((u2(1)/u1(1))^2) );
        % Outlet boundary condition
        u1(n)   = 2*u1(n-1) - u1(n-2);
        u2(n)   = 2*u2(n-1) - u2(n-2);
        u3(n)   = 2*u3(n-1) - u3(n-2);

        % Fluid properties at (k+1)th time step
        rho = u1./A;
        v   = u2./u1;
        T   = (g-1)*((u3./u1) - (g/2)*((u2./u1).^2));

        % Storing timewise variation of variables
        rho_t(:,:,k)   = rho;
        v_t(:,:,k)     = v;
        T_t(:,:,k)     = T;
        p_t(:,:,k)     = rho.*T;
        M_t(:,:,k)     = v./(T.^(0.5));
        mflow_t(:,:,k) = rho.*A.*v;


        % Checking for convergence
        if (k~=1)
           kn = k;              % Returning last timestep value
           error = max(abs(mflow_t(1,throat,k) - mflow_t(1,throat,k-1)));
           if (error < tol) 
           break
           end
        end

    end
end