clc;
clear;
close all;

%% VARIABLES
Hp = 0;          % [ft]
ISA = 0;             % [C or K]
speed_type = "mach"; % [Mach, KTAS, KEAS, KCAS]
speed_value = 0.7;   % [kn]


%% CALCULATIONS
atm = atmosphere_model(speed_type, speed_value, Hp, ISA);

disp(table(Hp, atm.T, atm.theta, atm.P, atm.delta, atm.rho, atm.sigma, atm.a, ...
    'VariableNames',{'Hp (ft)','T (K)','theta','P (N/m2)','delta','rho (kg/m3)','sigma','a (m/s)'}));

disp(table(atm.M, atm.KTAS, atm.KEAS, atm.KCAS, atm.q, atm.Rel, atm.mu, ...
    'VariableNames',{'Mach','KTAS','KEAS', 'KCAS','q (N/m2)','Re/l','mu (Ns/m2)'}));


%% PLOTTING
Hp_array = (1:500:160e3);

for i = 1:length(Hp_array)
    atm_array(i) = atmosphere_model(speed_type, speed_value, Hp_array(i), ISA);
end

figure
tiledlayout(2,4,'TileSpacing','Compact')
sgtitle( "Atmospheric Parameters")

nexttile
plot([atm_array.T],Hp_array)
grid minor
hold on
plot(atm.T,Hp, 'r.','MarkerSize', 15)
xlabel("Temperature")
ylabel("Altitude")

nexttile
plot([atm_array.P],Hp_array)
grid minor
hold on
plot(atm.P,Hp, 'r.','MarkerSize', 15)
xlabel("Pressure")
ylabel("Altitude")

nexttile
plot([atm_array.rho],Hp_array)
grid minor
hold on
plot(atm.rho,Hp, 'r.','MarkerSize', 15)
xlabel("Density")
ylabel("Altitude")

nexttile
plot([atm_array.a],Hp_array)
grid minor
hold on
plot(atm.a,Hp, 'r.','MarkerSize', 15)
xlabel("Speed of Sound")
ylabel("Altitude")

nexttile
plot([atm_array.M],Hp_array)
grid minor
hold on
plot(atm.M,Hp, 'r.','MarkerSize', 15)
xlabel("Mach Number")
ylabel("Altitude")

nexttile
plot([atm_array.KTAS],Hp_array)
grid minor
hold on
plot(atm.KTAS,Hp, 'r.','MarkerSize', 15)
xlabel("KTAS")
ylabel("Altitude")

nexttile
plot([atm_array.KEAS],Hp_array)
grid minor
hold on
plot(atm.KEAS,Hp, 'r.','MarkerSize', 15)
xlabel("KEAS")
ylabel("Altitude")

nexttile
plot([atm_array.KCAS],Hp_array)
grid minor
hold on
plot(atm.KCAS,Hp, 'r.','MarkerSize', 15)
xlabel("KCAS")
ylabel("Altitude")