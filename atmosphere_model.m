function atm = atmosphere_model(type, value, Hp, ISA)
%% INPUTS
% type: Velocity type, mach, KTAS, KEAS, KCAS
% value: Magnitude of velocity
% Hp: altitude [ft]
% ISA: Temperature difference [K]

%% OUTPUTS
% T: Temperature [K]
% P: Pressure [N/m2]

% function [AtmObject] = vel_model(type, value, Hp, ISA)
% fminbnd(@(hp) abs(vel_model('kcas',200,hp,0).KTAS - vel_model('mach',0.7,hp,0).KTAS),0,100e3)


%% ATMOSPHERE

    atm = struct();
    Hp = Hp * 0.3048;
    T_0 = 288.15;
    P_0 = 101325;
    rho_0 = 1.225;
    a_0 = 340.294 * 1.9438452;

    gamma = 1.4;
    R = 287.05287;

    if (0<=Hp) && (Hp<11e3)
        T = 288.15 + -6.5e-3 * Hp + ISA;
        P = (8.9619638 + -0.20216125e-3 * Hp) ^ (5.2558797);
        if ISA == 0
            rho = (1.048840 + -23.659414e-6 * Hp) ^ (4.2558797);
        else
            rho = P / (R * T);
        end

    elseif (11e3<=Hp) && (Hp<20e3)
        T = 216.65 + 0 * Hp + ISA;
        P = 128244.5 * exp(-0.15768852e-3 * Hp);
        if ISA == 0
            rho = 2.0621400 * exp(-0.15768852e-3 * Hp);
        else
            rho = P / (R * T);
        end

    elseif (20e3<=Hp) && (Hp<32e3)
        T = 196.65 + 1e-3 * Hp + ISA;
        P = (0.70551848 + 3.5876861e-6 * Hp) ^ (-34.163218);
        if ISA == 0
            rho = (0.9726309 + 4.94600e-6 * Hp) ^ (-35.163218);
        else
            rho = P / (R * T);
        end

    elseif (32e3<=Hp) && (Hp<47e3)
        T = 139.05 + 2.8e-3 * Hp + ISA;
        P = (0.34926867 + 7.0330980e-6 * Hp) ^ (-12.201149);
        if ISA == 0
            rho = (0.84392929 + 16.993902e-6 * Hp) ^ (-13.201149);
        else
            rho = P / (R * T);
        end

    elseif (47e3<=Hp) && (Hp<50e3)
        T = 270.65 + 0 * Hp + ISA;
        P = 41828.420 * exp(-0.12622656 * Hp);
        if ISA == 0
            rho = 0.53839563 * exp(-0.12622656 * Hp);
        else
            rho = P / (R * T);
        end

    else
        T = NaN;
        P = NaN;
        rho = NaN;
    end

    a = sqrt(gamma * R * T) * 1.9438452;

    theta = T / T_0;
    delta = P / P_0;
    sigma = rho / rho_0;
    mu = (1.458e-6 * T^(3/2)) / (T + 110.4);


    %% VELOCITY

    switch lower(type)
        %% MACH TO KTAS-KEAS-KCAS CONVERTIONS
        case 'mach'
            M = value;
            KTAS = value * a;
            KEAS = a_0 * value * sqrt(delta);
            q = 1.4 * P * (M^2) /2;

            if M <= 1
                P_p = P * (1+0.2*M^2)^3.5;
            elseif M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
            end

            p_k = (P_p - P) / P_0;

            if p_k <= 0.89293 && M <= 1
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k <= 0.89293 && M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
                p_k = (P_p - P) / P_0;
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k > 0.89293 && M > 1
                KCAS = a_0 * (.41726 + .7767*p_k - .0989 / p_k)^.5;
            end

        %% KTAS TO MACH-KEAS-KCAS CONVERTIONS
        case 'ktas'
            M = value / a;
            KTAS = value;
            KEAS = value  * sqrt(sigma);
            q = rho * KTAS / 2;

            if M <= 1
                P_p = P * (1+0.2*M^2)^3.5;
            elseif M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
            end

            p_k = (P_p - P) / P_0;

            if p_k <= 0.89293 && M <= 1
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k <= 0.89293 && M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
                p_k = (P_p - P) / P_0;
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k > 0.89293 && M > 1
                KCAS = a_0 * (.41726 + .7767*p_k - .0989 / p_k)^.5;
            end

        %% KEAS TO MACH-KTAS-KCAS CONVERTIONS
        case 'keas'
            M = value / (a_0 * sqrt(delta));
            KTAS = value / sqrt(sigma);
            KEAS = value;
            q = rho * KTAS / 2;

            if M <= 1
                P_p = P * (1+0.2*M^2)^3.5;
            elseif M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
            end

            p_k = (P_p - P) / P_0;

            if p_k <= 0.89293 && M <= 1
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k <= 0.89293 && M > 1
                P_p = P * ( (1.2*M^2)^3.5 * (1 + 7/6*(M^2-1))^-2.5 );
                p_k = (P_p - P) / P_0;
                KCAS = a_0 * (5*(p_k+1)^(2/7) - 5)^.5;
            elseif p_k > 0.89293 && M > 1
                KCAS = a_0 * (.41726 + .7767*p_k - .0989 / p_k)^.5;
            end

        %% KCAS TO MACH-KTAS-KEAS CONVERTIONS
        case 'kcas'
            M = sqrt(((((P_0 * ((((value / a_0) ^ 2 / 5) + 1) ^ (7/2) -1)) / P) + 1) ^ (2/7) - 1) / 0.2);
            
            %if (value / a_0) <= 1
            %    P_p = P_0 * (((1 + (0.2 * (value / a_0)^2)) ^ 3.5) - 1 + P);
            %elseif (value / a_0) > 1
            %    P_p = P_0 * ((((1.2 * (value / a_0)^2) ^ 3.5) * ((1 + (7/6) * ((value / a_0)^2 - 1)) ^ 2.5)) - 1 + P);
            %end
            %
            %p_k = P_p / P;
            %
            %if p_k <= 1.89293
            %    M = sqrt(5 * ((p_k ^ (2/7)) - 1));
            %elseif p_k > 1.89293
            %    M = sqrt((0.41726 + 0.7767 * (p_k - 1) - 0.0989 / (p_k - 1)));
            %end
             
            KTAS = M * a;
            KEAS = a_0 * M * sqrt(delta);
            KCAS = value;
            q = 1.4 * P * (M^2) /2;
    
        otherwise
            M = NaN;
            KTAS = NaN;
            KEAS = NaN;
            KCAS = NaN;
            q = NaN;
    end

    Rel = rho * a * value / mu;

    atm.T_0 = T_0;
    atm.P_0 = P_0;
    atm.rho_0 = rho_0;
    atm.a_0 = a_0;
    atm.T = T;
    atm.P = P;
    atm.rho = rho;
    atm.a = a;
    atm.gamma = gamma;
    atm.R = R;
    atm.theta = theta;
    atm.delta = delta;
    atm.sigma = sigma;
    atm.mu = mu;

    atm.M = M;
    atm.KTAS = KTAS;
    atm.KEAS = KEAS;
    atm.KCAS = KCAS;
    atm.q = q;
    atm.Rel = Rel;
    atm.mu = mu;
end