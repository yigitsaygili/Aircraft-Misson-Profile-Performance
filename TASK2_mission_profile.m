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

% MATRIX DISPLAY ---------------------------------------------------------
figure
tiledlayout(1,2)

nexttile
heatmap(alt_vector,mach_vector,thrust_grid, 'Colormap',autumn,'ColorbarVisible','off')
title("Thrust [N]")
xlabel("Altitude [ft]")
ylabel("Mach")

nexttile
heatmap(alt_vector,mach_vector,sfc_grid, 'Colormap',winter,'ColorbarVisible','off')
title("SFC [g/kN/s]")
xlabel("Altitude [ft]")
ylabel("Mach")


%% CALCULATIONS
cruise_alt_ft = 35000;
fprintf('Initial fuel : %.2f lbs\n', weights.max_fuel*2.2046226218488);
reserved_fuel = cruise_reserved_fuel(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid);

takeoff = takeoff_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid);
climb = climb_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, cruise_alt_ft, takeoff);
cruise = cruise_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, cruise_alt_ft, climb, reserved_fuel);
descent = descent_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, cruise_alt_ft, cruise);
landing = landing_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, descent);


%% RESULTS
total_distance = takeoff.total_distance/1000 + climb.total_distance + cruise.total_distance + descent.total_distance + landing.total_distance/1000;
total_time = (takeoff.total_time + climb.total_time + cruise.total_time*60 + descent.total_time + landing.total_time)/60;
altitude_path = [takeoff.altitude', climb.altitude, cruise.altitude', descent.altitude, landing.altitude'];

fprintf('Total flight distance: %.2f NM\n', total_distance*0.539956803);
fprintf('Total flight time    : %.2f h\n', total_time);

figure
plot(linspace(0,total_distance*0.539956803,length(altitude_path)),altitude_path,'r','LineWidth',1.5)
grid minor
title("Flight Trajectory")
xlabel("Horizontal Distance [NM]")
ylabel("Altitude [ft]")
pbaspect([4 1 1])
