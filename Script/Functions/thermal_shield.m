function [h_Tmax,Tmax] = thermal_shield(v,h,parT,T_w, q_rad_TS,q_cond,parT)

    %TO BE MODIFIED
    function [x,T,q_flux_in_Wcm2,dthabdt_cms] = Ablative_shield(N,parT)
        % Thermal protection system
        %
        % Inputs:
        % N: number of nodes in the solution (-)
        % Outputs:
        % x: position of nodes (m)
        % T: temperatures at nodes (K)
        % q_flux_in_Wcm^: heat flux to air (W/cm^2)
        % dthabdt_cms: rate of shield consumption (cm/s)

        th_ab = 0.036;              % ablation shield thickness (m)%From International Space Station (ISS) Soyuz Vehicle
                                    % Descent Module Evaluation of Thermal Protection
                                    % System (TPS) Penetration Characteristics
        th_s = 0.010;               % Steel thickness (m)
        k_s = 20;                   % thermal conductivity of steel (W/m-K) %MUST TO BE DETERMINATED
        q_flux = 100*100^2;         % heat flux (W/m^2)
        T_m = 755;                  % melting temperature ablative shield (K)
        DELTAi_fus_ab = 200e3;      % latent heat of fusion (J/kg)
        h_bar = 10;                 % heat transfer coefficient (W/m^2-K)
        T_inside = 320;             % internal air temperature (K)
        rho_ab = 1200;              % density (kg/m^3)
       
        % Setup nodes
        DELTAx = th_ab/(N-1);       %distance between adjacent nodes

        for i=1:N
            x(i,1) = th_ab*(i-1)/(N-1); %position of each node
        end

        % Resistances
        R_cond_s = th_s/(parT.A_ref*k_s);
        R_conv = 1/(parT.A_ref*h_bar);

        % Initial guess for temperature distribution
        for i =1:N
            Tg(i,1) = T_m + (T_inside - T_m)*(i-1)/(N-1); % Linear from melting to air
        end

        % Set up matrix A and vector b
        A = spalloc(N,N,3*N);         % Allocate space for sparse matrix.
                                        % This creates N-by-N all zero sparse
                                        % matrix with room to eventually hold 3N
                                        % nonzeros
        b = zeros(N,1);

        % Note: We have the equation T_1 = T_m. Here, in T_1 [1], the [1] indicates 
        % A(1,1) and  T_m indicates b(1)


        err = 0.999;                    % Initial value of error (K), must be larger
                                        % than tol tolerance for convergence
        tol = 0.01;
        while (err>tol)
            % Node 1
            A(1,1) = 1;
            b(1,1) = T_m;

            % Internal nodes
            for i=2:(N-1)
                A(i,i) = -k_ab((Tg(i) + Tg(i+1))/2)*A_c/DELTAx - k_ab((Tg(i) + Tg(i-1))/2)*A_c/DELTAx;
                A(i,i+1) = k_ab((Tg(i) + Tg(i+1))/2)*A_c/DELTAx;
                A(i,i-1) = k_ab((Tg(i) + Tg(i-1))/2)*A_c/DELTAx;
            end

            % Node N
            A(N,N) = -k_ab((Tg(N) + Tg(N-1))/2)*A_c/DELTAx - 1/(R_cond_s + R_conv);
            A(N,N-1) = k_ab((Tg(N) + Tg(N-1))/2)*A_c/DELTAx;
            b(N,1) = -T_inside/(R_cond_s + R_conv);


            % Solve matrix equation
            X = A\b;
            T = X;

            % Compute error
            err = sqrt(sum((T - Tg).^2)/N) 
            %Reset guess values used to setup A and b
            Tg = T;
        end


        % Calculation of required parameters

        q_flux_in = (T(N) - T_inside)/(R_cond_s + R_conv)/A_c;    % heat flux to air (W/m^2)
        q_flux_in_Wcm2 = q_flux_in/100^2;                           % heat flux to air (W/cm^2)
        dthabdt =(q_flux - q_flux_in)/(DELTAi_fus_ab*rho_ab);         % consumption rate (m/s)
        dthabdt_cms = dthabdt/100;
            end

end
