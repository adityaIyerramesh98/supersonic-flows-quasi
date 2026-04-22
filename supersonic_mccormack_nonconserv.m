function [rho_t,v_t,T_t,p_t,M_t,mflow_t,kn] = non_conservation_form(n,nt,C,gamma,tol)
x = linspace(0,3,n);              % Space Grid
dx  = x(2) - x(1);                % Space-step
A   = 1 + 2.2*(x - 1.5).^2;       % Area profile (Non-Dimensional)
% getting throat index
if (mod(n,2)==0),throat = n/2; % n is even
else ,throat = (n+1)/2;        % n is odd         
end

% Initial Conditions
rho = 1 - 0.3146*x;               % Density (Non-Dimensional)
T   = 1 - 0.2314*x;               % Temperature (Non-Dimensional)
v   = (0.1 + 1.09*x).*(T.^(0.5)); % Velocity (Non-Dimensional)

% Pre-allocating arrays for speed
[drhodt,dvdt,dTdt,drhodt_p,dvdt_p,dTdt_p] = deal(zeros(1,n-1));
[rho_t(:,:,1),v_t(:,:,1),T_t(:,:,1),...
    p_t(:,:,1),M_t(:,:,1),mflow_t(:,:,1)] = deal(zeros(1,n));
 
% Non-Conservation form of Governing Equ.
    for k = 1:nt
        dt  = min((C*dx)./((T.^0.5)+v));  % Selecting min value of dt 
                                          % across all grid points
        % Updating old values
        rho_old = rho;
        T_old = T;
        v_old = v;

        % Predictor Step
        for i = 2:n-1
            % Forward differences (at ith grid point)
            dvdx    = (v(i+1) - v(i))/dx;
            dlogAdx = (log(A(i+1)) - log(A(i)))/dx;
            drhodx  = (rho(i+1) - rho(i))/dx;
            dTdx    = (T(i+1) - T(i))/dx;
            
            % Time derivatives at (k)th time level
            % Continuity Equation
            drhodt(i) = -rho(i)*dvdx - rho(i)*v(i)*dlogAdx - v(i)*drhodx;
            % Momentum Equation
            dvdt(i)   = -v(i)*dvdx - (1/gamma)*(dTdx + (T(i)/rho(i))*drhodx);
            % Energy Equation
            dTdt(i)   = -v(i)*dTdx - (gamma-1)*T(i)*(dvdx + v(i)*dlogAdx);

            % Solution Update (Predicted values at (k+1)th time level)
            rho(i)  = rho_old(i) + drhodt(i)*dt;
            v(i)    = v_old(i) + dvdt(i)*dt;
            T(i)    = T_old(i) + dTdt(i)*dt;
        end

        % Corrector Step
        for j = 2:n-1
            % Backward differences (at jth grid point)
            dvdx    = (v(j) - v(j-1))/dx;
            dlogAdx = (log(A(j)) - log(A(j-1)))/dx;
            drhodx  = (rho(j) - rho(j-1))/dx;
            dTdx    = (T(j) - T(j-1))/dx;
            
            % Predicted time derivatives at (k+1)th time level
            % Continuity Equation
            drhodt_p(j) = -rho(j)*dvdx - rho(j)*v(j)*dlogAdx - v(j)*drhodx;
            % Momentum Equation
            dvdt_p(j)   = -v(j)*dvdx - (1/gamma)*(dTdx + (T(j)/rho(j))*drhodx);
            % Energy Equation
            dTdt_p(j)   = -v(j)*dTdx - (gamma-1)*T(j)*(dvdx + v(j)*dlogAdx);
        end

        % Average Time derivatives
        drhodt_av = 0.5*(drhodt + drhodt_p);
        dvdt_av   = 0.5*(dvdt + dvdt_p);
        dTdt_av   = 0.5*(dTdt + dTdt_p);
        
        % Final Values at time step k+1
        for m = 2:n-1
            rho(m) = rho_old(m) + drhodt_av(m)*dt;
            v(m)   = v_old(m) + dvdt_av(m)*dt;
            T(m)   = T_old(m) + dTdt_av(m)*dt;
        end
        
        % Boundary Conditions
        % Inlet boundary condition
        v(1)   = 2*v(2) - v(3);
        rho(1) = 1;
        T(1)   = 1;
        % Outlet boundary condition
        v(n)   = 2*v(n-1) - v(n-2);
        rho(n) = 2*rho(n-1) - rho(n-2);
        T(n)   = 2*T(n-1) - T(n-2);
        
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