clc
clear
close all


%% CONSTANTS
% WEIGHT DATA ------------------------------------------------------------
weights.MRW = 54700;
weights.MTOW = 54500;
weights.MLW = 45800;
weights.MZFW = 36500;
weights.empty = 32133;
weights.max_payload = 4366;
weights.payload_w_fuel = 694;
weights.max_fuel = 21871;

% AERODYNAMIC DATA -------------------------------------------------------
aero.Cd_0 = 0.0161;
aero.K = 0.0575;
aero.Cl_AoA_0 = 0.361;
aero.Cl_AoA_10 = 0.750;
aero.Cd_LG = 0.0150;
aero.Cd_AoA_10 = 0.0150;
aero.S = 92.5;
aero.FI_thrust = 0.8;
aero.GI_thrust = 0.6;
aero.rho = 1.225;
aero.g = 9.81;
aero.V = 80 * 0.5144;
aero.mu = 0.03;


%% THRUST - SFC DATA
% BREAKPOINTS ------------------------------------------------------------
alt_vector = [0, 18.5e3, 35e3, 41e3];
mach_vector = [0:0.1:0.7, 0.78, 0.85];

% THRUST MATRIX ----------------------------------------------------------
thrust_grid = [163960 96631 58184 44203;
               133286 86020 55744 44375;
               111052 76640 53347 44878;
               95840 68491 50995 44632;
               86233 61574 48686 43999;
               80814 55889 46421 42977;
               78166 51436 44199 41567;
               76872 48214 42021 39769;
               75859 46523 40310 38051;
               74368 45690 38836 36344];

% FUEL CONSUMPTION MATRIX ------------------------------------------------
sfc_grid = [14.5 14.3 14.1 13.9;
            15.3 14.8 14.4 14.1;
            15.8 15.3 14.9 14.6;
            16.2 15.8 15.3 15.0;
            16.7 16.2 15.8 15.5;
            17.2 16.7 16.2 15.9;
            17.7 17.2 16.7 16.4;
            18.2 17.7 17.2 16.8;
            18.7 18.1 17.6 17.3;
            19.1 18.5 18.0 17.6];


%% CALCULATIONS
% VARIABLES --------------------------------------------------------------
alt_ft = 0;
ISA = 15;
mass = weights.MTOW;
Weight_N = mass  * aero.g;
Lift_N = 0;

% FUEL AND WEIGHT CALCULATIONS -------------------------------------------
Wfuel_N(1) = weights.max_fuel * aero.g;
W1 = (weights.MTOW - weights.max_fuel) * aero.g;
Weight_N(1) = W1 + Wfuel_N(1);
mass_kg = Weight_N(1)/aero.g;

% INITITAL PARAMETERS ----------------------------------------------------
dt = 0.1;

v_kcas = 80;
t(1) = 0;
x_m(1) = 0;
V_mps(1) = 0;
i = 1;

% EULER INTEGRATION LOOP - 1 ---------------------------------------------
while V_mps<=v_kcas/1.9438452
    % THRUST CALCULATION -------------------------------------------------
    mach = atmosphere_model("kcas",V_mps(i)*1.9438452,alt_ft,ISA).M;
    Thrust_N = interpn(mach_vector, alt_vector, thrust_grid, mach, alt_ft);

    % FUEL CONSUMPTION CALCULATION ---------------------------------------
    sfc = interpn(mach_vector, alt_vector, sfc_grid, mach, alt_ft)/10e6;
    dW = -sfc * Thrust_N * dt * aero.g;
    
    Wfuel_N(i+1) = Wfuel_N(i) + dW;
    Weight_N(i+1) = W1 + Wfuel_N(i+1);
    mass_kg = Weight_N(i+1) / aero.g;
    
    % LIFT AND DRAG CALCULATION ------------------------------------------
    Lift_N = 0.5 * aero.rho * V_mps(i)^2 * aero.S * aero.Cl_AoA_0;
    Drag_N = 0.5 * aero.rho * V_mps(i)^2 * aero.S * (aero.Cd_LG + aero.Cd_0 + aero.K * aero.Cl_AoA_0^2);
    
    % CHANGE IN VELOCITY -------------------------------------------------
    dv = (Thrust_N - Drag_N - aero.mu * (Weight_N(i) - Lift_N)) / mass_kg;

    % ITERATIONS ---------------------------------------------------------
    t(i+1) = t(i) + dt; 
    V_mps(i+1) = V_mps(i) + dv*dt;
    x_m(i+1) = x_m(i) + V_mps(i)*dt + 0.5*dv*dt^2;
    
    i = i+1;
end

acc_t = t;
acc_x_m = x_m;
acc_V_mps = V_mps;

% EULER INTEGRATION LOOP - 2 ---------------------------------------------
while V_mps>=0
    % THRUST CALCULATION -------------------------------------------------
    mach = atmosphere_model("kcas",V_mps(i)*1.9438452,alt_ft,ISA).M;
    Thrust_N = interpn(mach_vector, alt_vector, thrust_grid, mach, alt_ft);

    % FUEL CONSUMPTION CALCULATION ---------------------------------------
    sfc = interpn(mach_vector, alt_vector, sfc_grid, mach, alt_ft)/10e6;
    dW = -sfc * Thrust_N * dt * aero.g;
    
    Wfuel_N(i+1) = Wfuel_N(i) + dW;
    Weight_N(i+1) = W1 + Wfuel_N(i+1);
    mass_kg = Weight_N(i+1) / aero.g;
    
    % LIFT AND DRAG CALCULATION ------------------------------------------
    Lift_N = 0.5 * aero.rho * V_mps(i)^2 * aero.S * aero.Cl_AoA_0;
    Drag_N = 0.5 * aero.rho * V_mps(i)^2 * aero.S * (aero.Cd_LG + aero.Cd_0 + aero.K * aero.Cl_AoA_0^2);
    
    % CHANGE IN VELOCITY -------------------------------------------------
    dv = (Thrust_N*0.06 - Drag_N - aero.mu*5 * (Weight_N(i) - Lift_N)) / mass_kg;

    % ITERATIONS ---------------------------------------------------------
    t(i+1) = t(i) + dt; 
    V_mps(i+1) = V_mps(i) + dv*dt;
    x_m(i+1) = x_m(i) + V_mps(i)*dt + 0.5*dv*dt^2;
    
    i = i+1;
end


%% RESULTS
total_distance = x_m(end);
total_time = t(end);
consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
reserved_fuel = Wfuel_N(end) / aero.g;

fprintf('Runway distance: %.2f m\n', total_distance);
fprintf('Total Time     : %.2f s\n', total_time);
fprintf('Fuel consumed  : %.2f kg\n', consumed_fuel);
fprintf('Fuel remaining : %.2f kg', reserved_fuel);

%% PLOTTING
figure
tiledlayout(1,2)
sgtitle('Aborted Takeoff Performance')

nexttile
plot(t, x_m, 'LineWidth',1.5)
hold on
plot(acc_t, acc_x_m, 'r', 'LineWidth',1.5)
hold off
grid minor
title("Distance - Time")
xlabel("Time [s]")
ylabel("Distance [m]")
legend("Decceleration", "Acceleration", 'Location','southeast')

nexttile
plot(t, V_mps*1.9438452, 'LineWidth',1.5)
hold on
plot(acc_t, acc_V_mps*1.9438452, 'r', 'LineWidth',1.5)
hold off
grid minor
title("Velocity - Time")
xlabel("Time [s]")
ylabel("Velocity [kn]")
legend("Decceleration", "Acceleration", 'Location','southeast')
