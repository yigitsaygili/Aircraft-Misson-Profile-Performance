function takeoff = takeoff_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid)


%% CALCULATIONS
% DEĞİŞKENLER ---------------------------------------------------------
previous_fuel_remaining = weights.max_fuel;
alt_ft = 0;
ISA = 15;
mass = weights.MTOW;
Weight_N = mass  * aero.g;
Lift_N = 0;

% YAKIT VE AĞIRLIK HESAPLAMALARI -----------------------------------------
Wfuel_N(1) = (previous_fuel_remaining) * aero.g;  % Yakıt ağırlığı [N]
W1 = (weights.MTOW - weights.max_fuel) * aero.g;  % Uçağın yakıtsız ağırlığı [N]
Weight_N(1) = W1 + Wfuel_N(1);  % Toplam ağırlık [N]
mass_kg = Weight_N(1)/aero.g;

% BAŞLANGIÇ PARAMETRELERİ ---------------------------------------------
dt = 0.1;
t(1) = 0;
x_m(1) = 0;
V_mps(1) = 0;

i = 1;
while Weight_N(i)-Lift_N>=0
    % İTKİ HESABI ------------------------------------------------------
    mach = atmosphere_model("ktas",V_mps(i)*1.9438452,alt_ft,ISA).M;
    Thrust_N = interpn(mach_vector, alt_vector, thrust_grid, mach, alt_ft);

    % YAKIT TÜKETİMİ HESAPLAMALARI ---------------------------------------
    sfc = interpn(mach_vector, alt_vector, sfc_grid, mach, alt_ft)/10e6;
    dW = -sfc * Thrust_N * dt * aero.g;
    
    Wfuel_N(i+1) = Wfuel_N(i) + dW;
    Weight_N(i+1) = W1 + Wfuel_N(i+1);
    mass_kg = Weight_N(i+1) / aero.g;
    
    % TAŞIMA - SÜRÜKLEME HESABI ---------------------------------------
    Lift_N_AoA_0 = 0.5 * aero.rho * V_mps(i)^2 * aero.S * aero.Cl_AoA_0;
    Lift_N_AoA_10 = 0.5 * aero.rho * V_mps(i)^2 * aero.S * aero.Cl_AoA_10;
    
    if Weight_N(i)-Lift_N_AoA_10<=0
        Lift_N = Lift_N_AoA_10;
    else
        Lift_N = Lift_N_AoA_0;
    end

    Drag_N = 0.5 * aero.rho * V_mps(i)^2 * aero.S * (aero.Cd_LG + aero.Cd_0 + aero.K * aero.Cl_AoA_0^2);
  
    % HIZ DEĞİŞİMİ HESABI ---------------------------------------------
    dv = (Thrust_N - Drag_N - aero.mu * (Weight_N(i) - Lift_N)) / mass_kg;

    % İTERASYONLAR ----------------------------------------------------
    t(i+1) = t(i) + dt; 
    V_mps(i+1) = V_mps(i) + dv*dt;
    x_m(i+1) = x_m(i) + V_mps(i)*dt + 0.5*dv*dt^2;
    
    i = i+1;
end


%% SONUÇLAR
total_distance = x_m(end);
total_time = t(end) / 60;
consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
reserved_fuel = Wfuel_N(end) / aero.g;

fprintf("\nTAKE-OFF --------------------\n");
fprintf('Takeoff distance: %.2f ft\n', total_distance*3.280839895);
fprintf('Takeoff time    : %.2f min\n', total_time);
fprintf('Fuel consumed   : %.2f lbs\n', consumed_fuel*2.2046226218488);
fprintf('Fuel remaining  : %.2f lbs\n', reserved_fuel*2.2046226218488);

%% PLOTTING
figure%('Name','Takeoff-1')
tiledlayout(1,2)
sgtitle("Takeoff Performance")

nexttile
plot(t, x_m)
grid minor
title("Distance - Time")
xlabel("Time [s]")
ylabel("Distance [m]")

nexttile
plot(t, V_mps)
grid minor
title("Velocity - Time")
xlabel("Time [s]")
ylabel("Velocity [kn]")

%% FUNCTION OUTPUT
takeoff.total_distance = x_m(end);
takeoff.total_time = t(end) / 60;
takeoff.consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
takeoff.reserved_fuel = Wfuel_N(end) / aero.g;
takeoff.final_velocity = V_mps(end);
takeoff.altitude = zeros(length(t),1);

end