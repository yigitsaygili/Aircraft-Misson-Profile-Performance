function cruise = cruise_performance(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid, cruise_alt_ft, climb, reserved_fuel)


%% INITIAL CONDITIONS
previous_fuel_remaining = climb.reserved_fuel;
ISA = 0;                       % ISA deviasyonu [K]
alt_ft = cruise_alt_ft;                % Seyir irtifası [ft]
dt = 1;                        % Zaman adımı [s]
mach = 0.78;                   % Sabit hız [KCAS]

% HIZ HESAPLAMALARI ------------------------------------------------------
tas_mps = atmosphere_model("mach", mach, alt_ft, ISA).KTAS / 1.9438452; % KTAS'ı m/s'ye çevir

% YAKIT VE AĞIRLIK HESAPLAMALARI -----------------------------------------
Wfuel_reserve_N = (reserved_fuel+ 10.767992 + 14.497105) * aero.g;

Wfuel_N(1) = (previous_fuel_remaining) * aero.g;  % Yakıt ağırlığı [N]
W1 = (weights.MTOW - weights.max_fuel) * aero.g;  % Uçağın yakıtsız ağırlığı [N]
Weight_N(1) = W1 + Wfuel_N(1);  % Toplam ağırlık [N]

% İLK DEĞERLER -----------------------------------------------------------
Range_km(1) = 0;  % Başlangıçta kat edilen mesafe [km]
t(1) = 0;  % Başlangıç zamanı [s]
i = 1;  % İndeks


%% CALCULATIONS - EULER
while Wfuel_N(i)>=Wfuel_reserve_N
    % YOĞUNLUK HESAPLAMASI -----------------------------------------------
    atmosphere_object = atmosphere_model("ktas", tas_mps * 1.9438452, alt_ft, ISA);
    rho = atmosphere_object.rho;

    % TAŞIMA VE SÜRÜKLEME HESAPLAMALARI --------------------------------------
    cl = (2 * Weight_N(i)) / (rho * tas_mps^2 * aero.S);
    Drag_N = 0.5 * rho * tas_mps^2 * aero.S * (aero.Cd_0 + aero.K * cl^2);

    % YAKIT TÜKETİMİ HESAPLAMALARI ---------------------------------------
    thrust_required = Drag_N;  % Sürüklenmeyi dengelemek için gereken itki
    sfc = interpn(mach_vector, alt_vector, sfc_grid, mach, alt_ft)/1e6; %[kg/Ns]
    dW = -sfc * thrust_required * dt * aero.g;  % Ağırlık azalması [N]
    
    Wfuel_N(i+1) = Wfuel_N(i) + dW;
    Weight_N(i+1) = W1 + Wfuel_N(i+1);
    mass_kg = Weight_N(i+1) / aero.g;

    % MESAFE HESAPLAMALARI -------------------------------------------------
    Range_km(i+1) = Range_km(i) + tas_mps * dt / 1000;  % Kat edilen mesafe [km]

    % ZAMAN GÜNCELLEMESİ ---------------------------------------------------
    t(i+1) = t(i) + dt;
    i = i + 1;
end


%% RESULTS
total_distance = Range_km(end);  % Toplam kat edilen mesafe [km]
total_time = t(end) / 3600;     % Toplam uçuş süresi [saat]
consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;    % Kalan yakıt miktarı [kg]
reserved_fuel = Wfuel_N(end) / aero.g;    % Kalan yakıt miktarı [kg]

fprintf("\nCRUISE --------------------\n");
fprintf('Cruise range   : %.2f NM\n', total_distance*0.539956803);
fprintf('Cruise time    : %.2f h\n', total_time);
fprintf('Fuel consumed  : %.2f lbs\n', consumed_fuel*2.2046226218488);
fprintf('Fuel remaining : %.2f lbs\n', reserved_fuel*2.2046226218488);
fprintf('Cruise velocity: %.2f mach\n', mach);
fprintf('Cruise alttude : %.2f ft\n', alt_ft);


%% PLOTTING
figure%('Name','Cruise-1')
tiledlayout(2,2)
sgtitle("Cruise Performance")

nexttile
plot(t/3600, Range_km,"k")
grid minor
title("Distance - Time")
ylabel("Distance [km]")
xlabel("Time [hours]")

nexttile
plot(t/3600, Weight_N/aero.g,"k")
grid minor
title("Aircraft Weight - Time")
ylabel("Weight [kg]")
xlabel("Time [hours]")

nexttile
plot(Range_km, Weight_N/aero.g,"k")
grid minor
title("Aircraft Weight - Distance")
ylabel("Weight [kg]")
xlabel("Distance [km]")

nexttile
plot(t/3600, Wfuel_N/aero.g,"k")
grid minor
title("Fuel Weight - Time")
ylabel("Fuel Weight [kg]")
xlabel("Time [hours]")


figure%('Name','Cruise-2')
tiledlayout(2,1)
sgtitle("Cruise Performance")

nexttile
plot(t/3600, alt_ft*(ones(length(t),1)))
grid minor
title("Altitude - Time")
ylabel("Altitıde [ft]")
xlabel("Time [hours]")

nexttile
plot(t/3600, tas_mps*(ones(length(t),1))*1.9438452)
grid minor
title("KTAS - Time")
ylabel("KTAS [kn]")
xlabel("Time [hours]")


%% FUNCTION OUTPUT
cruise.total_distance = Range_km(end);
cruise.total_time = t(end) / 3600;
cruise.consumed_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
cruise.reserved_fuel = Wfuel_N(end) / aero.g;
cruise.final_velocity = tas_mps(end);
cruise.altitude = alt_ft*(ones(length(t),1));

end