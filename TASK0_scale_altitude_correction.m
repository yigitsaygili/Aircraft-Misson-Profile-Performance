%clc;
clear;
close all

%% VARIABLES
ISA = 0;        % [C or K]
Hp_array = [0:2e3:60e3 65e3 : 5e3: 100e3];
Vc_array = (100:10:1000);

figure
axes('XAxisLocation','top')
hold on

%% CONSTANT ALTITUDE CURVES [ESDU 69026 TABLE 3A]
for j = 1:length(Hp_array)
    EAS_array = Vc_array*nan;
    CAS_array = Vc_array*nan;

    for i = 1:length(Vc_array)
        if (Hp_array(j) < 18e3 && Vc_array(i) < 200 ) || ...
           (Hp_array(j) < 40e3 && Vc_array(i) < 175 ) || ...
           (Hp_array(j) < 50e3 && Vc_array(i) < 150 ) || ...
           (Hp_array(j) < 60e3 && Vc_array(i) < 125 )
            continue
        end
        [~,~, EAS_array(i), CAS_array(i),~,~,~] = vel_model('kcas', Vc_array(i), Hp_array(j), ISA);
    end

    dV = EAS_array - CAS_array;
    plot(CAS_array, dV, 'k')
end


%% CONSTANT MACH NUMBER CURVES [ESDU 69026 TABLE 3A]

clear EAS_array CAS_array
M_array = (0.5 : .1 : 2.5);

for j = 1:length(M_array)
    Vc_arrayOut = M_array * nan;

    for i = 1:length(Hp_array)
        [~,~, atm.EAS_array(i), Vc_arrayOut(i),~,~,~] = vel_model('Mach', M_array(j), Hp_array(i), ISA);
    end

    dV = atm.EAS_array - Vc_arrayOut;
    plot(Vc_arrayOut, dV, 'k')
end

%% PLOT REFINEMENT
ylim([-55 0])
xlim([100 1000])
xlabel("Vc")
ylabel("Ve - Vc")
title("Scale Altitude Correction")
grid on
grid minor