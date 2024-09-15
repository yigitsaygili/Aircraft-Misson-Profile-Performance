function climb = climb_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, cruise_alt_ft, takeoff)


%% BAŞLANGIÇ HESAPLAMALARI
previous_fuel_remaining = takeoff.reserved_fuel;
ISA = 0;                       % ISA farkı [K]
alt_ft = 0;                % Başlangıç irtifası [ft]
v_kcas = 250;                  % Sabit hız [KCAS]
dt = 1;                        % Zaman adımı [s]
dH_ft = 500;                   % İrtifa adımı [ft] dV/dH hesabı için geçici parametre

% YAKIT VE AĞIRLIK HESAPLAMALARI -----------------------------------------
Wfuel_N(1) = (previous_fuel_remaining) * aero.g;  % Yakıt ağırlığı [N]
W1 = (weights.MTOW - weights.max_fuel) * aero.g;  % Uçağın yakıtsız ağırlığı [N]
Weight_N(1) = W1 + Wfuel_N(1);  % Toplam ağırlık [N]

% İLK DEĞERLER -----------------------------------------------------------
altitude_ft(1) = alt_ft;
distance_km(1) = 0;

tas_mps(1) = atmosphere_model("kcas", v_kcas, altitude_ft(1), ISA).KTAS / 1.9438452; % KTAS'ı m/s'ye çevir
theta(1) = 0;
roc_mps(1) = 0;
climb_mach(1) = atmosphere_model("ktas", tas_mps(1), altitude_ft(1), ISA).M;

t(1) = 0;
i = 1;


%% ALÇALMA HESAPLAMALARI (EULER YÖNTEMİ)
while altitude_ft(i) < cruise_alt_ft
    % YOĞUNLUK HESAPLAMASI -----------------------------------------------
    atm_object = atmosphere_model("ktas", tas_mps(i) * 1.9438452, altitude_ft(i), ISA);
    rho = atm_object.rho;

    % İTKİ HESAPLAMALARI -------------------------------------------------
    mach = atm_object.M; % V: kn
    Thrust_N = interpn(mach_vector, alt_vector, thrust_grid, mach, altitude_ft(i));

    % TAŞIMA VE SÜRÜKLEME HESAPLAMALARI --------------------------------------
    cl = (2 * Weight_N(i) * cos(theta(i))) / (rho * tas_mps(i)^2 * aero.S);
    Drag_N = 0.5 * rho * tas_mps(i)^2 * aero.S * (aero.Cd_0 + aero.K * cl^2);

    % YAKIT TÜKETİMİ HESAPLAMALARI ---------------------------------------
    sfc = interpn(mach_vector, alt_vector, sfc_grid, mach, altitude_ft(i))/10e6;
    dW = -sfc * Thrust_N * dt * aero.g;
    
    Wfuel_N(i+1) = Wfuel_N(i) + dW;
    Weight_N(i+1) = W1 + Wfuel_N(i+1);
    mass_kg = Weight_N(i+1) / aero.g;

    % RATE OF DESCENT ----------------------------------------------------
    atm_new_object = atmosphere_model("ktas", tas_mps(i)*1.9438452, altitude_ft(i)+dH_ft, ISA);
    dvdh =  ((atm_new_object.KCAS  -  atm_object.KCAS)* 0.514444) / (dH_ft / 0.3048);
    roc_mps(i+1) = (Thrust_N - Drag_N) / (mass_kg*dvdh + Weight_N(i)/tas_mps(i));
    theta(i+1) = asin(roc_mps(i+1)/tas_mps(i));

    % İRTİFA VE MESAFE HESAPLAMALARI --------------------------------------
    altitude_ft(i+1) = altitude_ft(i) + roc_mps(i) * dt / 0.3048;  % Alçalmış irtifa [ft]
    distance_km(i+1) = distance_km(i) + tas_mps(i) * dt / 1000;  % Kat edilen mesafe [km]
    tas_mps(i+1) = atmosphere_model("kcas", v_kcas, altitude_ft(i+1), ISA).KTAS / 1.9438452; % KTAS'ı m/s'ye çevir
    climb_mach(i+1) = atmosphere_model("ktas", tas_mps(i+1), altitude_ft(i+1), ISA).M;

    % ZAMAN GÜNCELLEMESİ ---------------------------------------------------
    t(i+1) = t(i) + dt;
    i = i + 1;
end


%% RESULTS
total_distance = distance_km(end);  % Toplam kat edilen mesafe [km]
total_time = t(end) / 60;  % Toplam alçalma süresi [saat]
consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;    % Kalan yakıt miktarı [kg]
reserved_fuel = Wfuel_N(end) / aero.g;    % Kalan yakıt miktarı [kg]

fprintf("\nCLIMB --------------------\n");
fprintf('Climb distance: %.2f NM\n', total_distance*0.539956803);
fprintf('Climb time    : %.2f min\n', total_time);
fprintf('Fuel consumed : %.2f lbs\n', consumed_fuel*2.2046226218488);
fprintf('Fuel remaining: %.2f lbs\n', reserved_fuel*2.2046226218488);
fprintf('Climb velocity: %.2f kn\n', v_kcas);
fprintf('Final altitude: %.2f ft\n', cruise_alt_ft);


%% PLOTTING
figure%('Name','Climb-1')
tiledlayout(2,2)
sgtitle("Climb Performance")

nexttile
plot(t/60, altitude_ft,"k")
grid minor
title("Altitude - Time")
ylabel("Altitude [ft]")
xlabel("Time [min]")

nexttile
plot(t/60, Weight_N/aero.g,"k")
grid minor
title("Aircraft Weight - Time")
ylabel("Weight [kg]")
xlabel("Time [min]")

nexttile
plot(t/60, roc_mps,"k")
grid minor
title("Rate of Descent - Time")
ylabel("ROD [m/s]")
xlabel("Time [min]")

nexttile
plot(t/60, theta,"k")
grid minor
title("Descent Angle - Time")
ylabel("Theta [rad]")
xlabel("Time [min]")


figure%('Name','Climb-2')
tiledlayout(1,2)
sgtitle("Climb Performance")

nexttile
plot(t/60, tas_mps*1.9438452)
grid minor
title("KTAS - Time")
ylabel("KTAS [kn]")
xlabel("Time [min]")

nexttile
plot(tas_mps*1.9438452, altitude_ft)
grid minor
title("Altitude - KTAS")
ylabel("Altitude [ft]")
xlabel("KTAS [kn]")




%% FUNCTION OUTPUT
climb.total_distance = distance_km(end);
climb.total_time = t(end) / 60;
climb.consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
climb.reserved_fuel = Wfuel_N(end) / aero.g;
climb.final_velocity = tas_mps(end);
climb.altitude = altitude_ft;

end