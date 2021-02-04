function plotAtmosphereV2(H_layer, lambda, Tb, Pb, R, g0, Mm)
%--------------------------------------------------------------------------
% Inputs:
%   - H_layer   Atmospheric layers altitude                         [m]
%   - lambda    Atmospheric layers thermal gradients                [K/m]
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - Pb        Vector of base pressures, 1 x length(H_layer)       [Pa]
%   - R         Universal gas constant                              [J/(kg K)]
%   - g0        Acceleration at planet's surface                    [m/s^2]
%   - Mm        Planet's air molar mass                             [kg/mol]
%--------------------------------------------------------------------------
% Outputs:
%       - Temperature vs. Altitude
%       - Pressure vs. Altitude (up until 20 km)
%       - Pressure vs. Altitude (from 20 km)
%       - Density vs. Altitude
%--------------------------------------------------------------------------

% Plot Temperature vs. Altitude
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Temperatura vs. Altitud}");
for i = 1:length(H_layer)-1
    plot([Tb(i) Tb(i+1)], [H_layer(i) H_layer(i+1)]/1e3, 'r');
end
xlabel("Temperatura $\left( \mathrm{K} \right)$");
ylabel("Altitud $\left( \mathrm{km} \right)$");
yticks(10*[0:2:12]);
grid on;
grid minor;
box on;
set(gcf, 'units', 'centimeters', 'position', [32,1,18,15]);
hold off;

% Plot Pressure vs. Altitude (until 20 km)
DeltaH = 1000;
H_plot = 0:DeltaH:max(H_layer);
P_plot = zeros(1, length(H_plot));
index_20 = 20e3/DeltaH + 1;
for i = 1:length(H_plot)
    P_plot(i) = getPressureV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H_plot(i));
    P_plot(i) = abs(P_plot(i));
    %fprintf("%15d %15.2f %15.2f\n", i, H_plot(i), P_plot(i));
end

figure(2);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Pressi\'o vs. Altitud}");
plot(P_plot(1:index_20)/1e2, H_plot(1:index_20)/1e3, 'b');
xlabel("Pressi\'o $\left( \mathrm{kPa} \right)$");
ylabel("Altitud $\left( \mathrm{km} \right)$");
yticks([0:2:20]);
grid on;
grid minor;
box on;
set(gcf, 'units', 'centimeters', 'position', [32,1,18,15]);
hold off;

% Plot Pressure vs. Altitude (from 20 kg)
figure(3);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Pressi\'o vs. Altitud}");
plot(P_plot(index_20+1:end)/1e2, H_plot(index_20+1:end)/1e3, 'b');
xlabel("Pressi\'o $\left( \mathrm{kPa} \right)$");
ylabel("Altitud $\left( \mathrm{km} \right)$");
grid on;
grid minor;
box on;
set(gcf, 'units', 'centimeters', 'position', [32,1,18,15]);
hold off;

% Plot Density vs. Altitude
rho_plot = zeros(1, length(H_plot));
for i = 1:length(H_plot)
    rho_plot(i) = getDensityV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H_plot(i));
    rho_plot(i) = abs(rho_plot(i));
    %fprintf("%15d %15.2f %15.2f\n", i, H_plot(i), P_plot(i));
end

figure(4);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Densitat vs. Altitud}");
plot(rho_plot, H_plot/1e3, 'b');
xlabel("Densitat $\left( \mathrm{kg} / \mathrm{m}^3 \right)$");
ylabel("Altitud $\left( \mathrm{km} \right)$");
grid on;
grid minor;
box on;
set(gcf, 'units', 'centimeters', 'position', [32,1,18,15]);
hold off;
hold off;

end