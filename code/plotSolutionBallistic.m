function plotSolutionBallistic(beta, H, V, gamma, t, r, a, P_dynamic, Ma, Re, ...
    P_stagnation, h_stagnation, HT_stagnation, E_dynamic, printTitle)
%--------------------------------------------------------------------------
% Inputs: 
%       - beta          Ballistic coefficient vectors       [Pa]
%       - H             Heights matrix                      [m]
%       - V             Velocities matrix                   [m/s]
%       - t             Times matrix                        [s]
%       - r             Ranges matrix                       [m]
%       - a             Acceleration matrix                 [g]
%       - P_dynamic     Dynamic pressure matrix             [Pa]
%       - Ma            Mach number matrix                  [adim]
%       - Re            Reynolds number matrix              [1/m]
%       - P_stagnation  Stagnation pressure matrix          [Pa]
%       - h_stagnation  Stagnation enthalpy matrix          [J/kg]
%       - HT_stagnation Stagnation heat transfer matrix     [W/m^2]
%       - E_dynamic     Dynamic energy matrix               [Pa·m/s]
%       - num_fig       Number of previous figure
%       - printTitle    Boolean argument to decide if print title or not
% Note: every solution matrix consist of n rows, one for each gamma0. In
% every row, the solution corresponding to the magnitude is stored
%--------------------------------------------------------------------------
% Outputs: 
%       - Altitude vs. Velocity plot
%       - Altitude vs. Acceleration plot
%       - Altitude vs. Dynamic pressure plot
%       - Altitude vs. Mach number plot
%       - Altitude vs. Reynolds number plot
%       - Altitude vs. Stagnation Point Pressure plot
%       - Altitude vs. Stagnation Point Enthalpy plot
%       - Altitude vs. Entry time plot
%       - Altitude vs. Range plot 
%       - Altitude vs. Dynamic energy plot
%--------------------------------------------------------------------------

% Colors
color = ['b', 'r', 'g', 'm', 'c', 'y', 'k'];
beta = beta/(unitsratio('feet', 'meter')^2*4.4482216152605);
legend_str = {sprintf('$\\beta = %.0f \\ %s$', beta(1), '\mathrm{lbf}/\mathrm{ft}^2')};
for i = 2:length(beta)
     legend_str(end+1) = {sprintf('$\\beta = %.0f \\ %s$', beta(i), '\mathrm{lbf}/\mathrm{ft}^2')}; %#ok<AGROW>
end

% Plot size and location
posX = 0;
posY = 1;
sizeX = 18;
sizeY = 15;

% Altitude vs. Velocity
figure('Name', 'Ballistic Reentry - Velocity vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Velocitat}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(V(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Velocitat $\left( \mathrm{m} / \mathrm{s} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northwest');
hold off;

% Deceleration vs. Altitude
figure('Name', 'Ballistic Reentry - Deceleration vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Desacceleraci\'o}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(a(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Desacceleraci\'o $\left( g_0 \right)$");
xlim([0 40]);
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Dynamic Pressure vs. Altitude
figure('Name', 'Ballistic Reentry - Dynamic Pressure vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Pressi\'o din\`amica}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(P_dynamic(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Pressi\'o din\`amica $\left( \mathrm{Pa} \right)$");
%ylim([0 7660841.41]);
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Mach number vs. Altitude
figure('Name', 'Ballistic Reentry - Mach Number vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{N\'umero de Mach}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(Ma(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("N\'umero de Mach");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northwest');
hold off;

% Reynolds number vs. Altitude
figure('Name', 'Ballistic Reentry - Reynolds Number vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{N\'umero de Reynolds}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(Re(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("N\'umero de Reynolds $\left( \mathrm{m}^{-1} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Stagnation Point Pressure vs. Altitude
figure('Name', 'Ballistic Reentry - Stagnation Point Pressure vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Pressi\'o d'estancament}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(P_stagnation(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Pressi\'o d'estancament $\left( \mathrm{Pa} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Stagnation Point Enthalpy vs. Altitude
figure('Name', 'Ballistic Reentry - Stagnation Point Enthalpy vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Entalpia d'estancament}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(h_stagnation(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Entalpia d'estancament $\left( \mathrm{J} / \mathrm{kg} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northwest');
hold off;

% Stagnation Point Heat Transfer vs. Altitude
figure('Name', 'Ballistic Reentry - Stagnation Point Heat Transfer vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Transfer\`encia de calor d'estancament}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(HT_stagnation(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Transfer\`encia de calor d'estancament $\left( \mathrm{W} / \mathrm{m^2} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northwest');
hold off;

% Time vs. Altitude
figure('Name', 'Ballistic Reentry - Time vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Temps}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(t(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Temps $\left( \mathrm{s} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Range vs. Altitude
figure('Name', 'Ballistic Reentry - Range vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Abast}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(r(i,:)/1e3, H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Abast $\left( \mathrm{km} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Dynamic Energy vs. Altitude
figure('Name', 'Ballistic Reentry - Dynamic Energy vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Energia din\`amica}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(E_dynamic(i,:), H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Energia din\`amica $\left( \mathrm{Pa \cdot m} / \mathrm{s} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;



% Flight path angle vs. Altitude
figure('Name', 'Lifting Reentry - Flight path angle vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Angle de descens vs. Altitud}");
set(groot,'defaultAxesTickLabelInterpreter','latex');
end
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(gamma(i,:)*180/pi, H/1e3, color(i));
end
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Angle de descens $\left( ^\circ \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northeast');
hold off;

% Flight path angle vs. Time
figure('Name', 'Lifting Reentry - Flight path angle vs. Time');
hold on;
if printTitle == 1
    title("\textbf{Angle de descens vs. Temps}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(beta)
    plot(t(i,:), gamma(i,:)*180/pi, color(i));
end
xlabel("Temps $\left( \mathrm{s} \right)$");
ylabel("Angle de descens $\left( ^\circ \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
legend(legend_str, 'location', 'northwest');
hold off;

end

