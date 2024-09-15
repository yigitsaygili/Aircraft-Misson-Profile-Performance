function reserved_fuel = cruise_reserved_fuel(weights, aero, alt_vector, mach_vector, thrust_grid, sfc_grid)


%% INITIAL CONDITIONS
previous_fuel_remaining = weights.max_fuel;
ISA = 0;                       % ISA deviasyonu [K]
alt_ft = 35000;                % Seyir irtifası [ft]
dt = 1;                        % Zaman adımı [s]
mach = 0.78;                   % Sabit hız [KCAS]

% HIZ HESAPLAMALARI ------------------------------------------------------
tas_mps = atmosphere_model("mach", mach, alt_ft, ISA).KTAS / 1.9438452; % KTAS'ı m/s'ye çevir

% YAKIT VE AĞIRLIK HESAPLAMALARI -----------------------------------------
Wfuel_N(1) = (previous_fuel_remaining) * aero.g;  % Yakıt ağırlığı [N]
W1 = (weights.MTOW - weights.max_fuel) * aero.g;  % Uçağın yakıtsız ağırlığı [N]
Weight_N(1) = W1 + Wfuel_N(1);  % Toplam ağırlık [N]

% İLK DEĞERLER -----------------------------------------------------------
Range_km(1) = 0;  % Başlangıçta kat edilen mesafe [km]
t(1) = 0;  % Başlangıç zamanı [s]
i = 1;  % İndeks


%% CALCULATIONS - EULER
while t(i)<=30*60
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


%% FUNCTION OUTPUT
reserved_fuel = (Wfuel_N(1)-Wfuel_N(end)) / aero.g;
fprintf('Reserved fuel: %.2f lbs\n', reserved_fuel*2.2046226218488);

end