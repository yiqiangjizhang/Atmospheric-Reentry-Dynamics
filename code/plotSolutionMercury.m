function plotSolutionMercury(H, V, t, r, a, Ma, printTitle)
%--------------------------------------------------------------------------
% Inputs: 
%       - H             Heights matrix                      [m]
%       - V             Velocities matrix                   [m/s]
%       - t             Times matrix                        [s]
%       - r             Ranges matrix                       [m]
%       - a             Acceleration matrix                 [g0]
%       - num_fig       Number of previous figure
%--------------------------------------------------------------------------
% Outputs: 
%       - Altitude vs. Velocity plot
%       - Altitude vs. Entry time plot
%       - Altitude vs. Range plot 
%       - Altitude vs. Acceleration plot
%       - Altitude vs. Mach number plot
%--------------------------------------------------------------------------

% Plot size and location
posX = 16;
posY = 1;
sizeX = 18;
sizeY = 15;

% Altitude vs. Velocity
figure('Name', 'Mercury Capsule Reentry - Velocity vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Velocitat}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(V, H/1e3, 'b');
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Velocitat $\left( \mathrm{m} / \mathrm{s} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
hold off;

% Time vs. Altitude
figure('Name', 'Mercury Capsule Reentry - Time vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Temps}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(t, H/1e3, 'b');
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Temps $\left( \mathrm{s} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
hold off;

% Range vs. Altitude
figure('Name', 'Mercury Capsule Reentry - Range vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Abast}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(r/1e3, H/1e3, 'b');
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Abast $\left( \mathrm{km} \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
hold off;

% Altitude vs. Acceleration
figure('Name', 'Mercury Capsule Reentry - Deceleration vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{Desacceleraci\'o}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(a, H/1e3, 'b');
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("Desacceleraci\'o $\left( g_0 \right)$");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
hold off;

% Mach number vs. Altitude
figure('Name', 'Mercury Capsule Reentry - Mach Number vs. Altitude');
hold on;
if printTitle == 1
    title("\textbf{N\'umero de Mach}");
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(Ma, H/1e3, 'b');
ylabel("Altitud $\left( \mathrm{km} \right)$");
xlabel("N\'umero de Mach");
grid on;
grid minor;
set(gcf, 'units', 'centimeters', 'position', [posX, posY, sizeX, sizeY]);
box on;
hold off;

end

