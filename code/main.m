%% ATMOSPHERIC RE-ENTRY

clear all;
close all;
clc;

%% 1. Definition of Constants, Parameters and Variables

% 1.1. CONSTANTS
    
    G = 6.67408e-11;            % Gravitational Constant                [m^3/(kg*s^2)]
    R = 8.31432;                % Universal Constant for Ideal Gases    [J/mole*K]
    
    % 1.1.1. Earth
    R_Earth = 6371.0e3;         % Radius of Earth                       [m]
    M_Earth = 5.9724e24;        % Earth's Mass                          [kg]
    g0_Earth = 9.80665;         % Acceleration at Earth's surface       [m/s^2]
    T0_Earth = 288.15;          % US Standard Sea Level Temperature     [K]
    P0_Earth = 101325;          % Pressure at Sea Level                 [Pa]
    Mm_Earth = 28.9644*10^-3;   % Molecular Mass                        [kg*mole^-1]    
    % Earth's atmospheric layers altitude           [m]
    H_layer_Earth = 1e3*[0 11 20 32 47 52 61 69 79 90 100 110 117.776];    
    % Earth's atmopheric layers thermal gradient    [K/m]    
    lambda_layer_Earth = 1e-3*[-6.5 0 1 2.8 0 -2 -4 -3 0 2 4.36 16.4596 0];
    gamma_gas_Earth = 1.4;      % Earth's air specific heats relation   [adim]
    
    % 1.1.2. Mars
    R_Mars = 3389.5e3;          % Mars Radius                           [m]
    M_Mars = 0.64171e24;        % Mars' Mass                            [kg]
    g0_Mars = 3.71;             % Acceleration at Mars' surface         [m/s^2]
    T0_Mars = 242.15;           % US Standard Sea Level Temperature     [K]
    P0_Mars = 636;              % Pressure at Sea Level                 [Pa]
    Mm_Mars = 28.9644*10^-3;    % Molecular Mass                        [kg*mole^-1]    
    % Mars' atmospheric layers altitude             [m]
    H_layer_Mars = 1e3*[0 7 117.776];   
    % Mars' atmopheric layers thermal gradient      [K/m]    
    lambda_layer_Mars = 1e-3*[-0.998 -2.22];
    gamma_gas_Mars = 1.28;      % Mars' air specific heats relation   [adim]
    
% 1. 3. VARIABLES
syms V;                     % Vehicle velocity                          [m/s]
syms gamma;                 % Vehicle flight path angle                 [rad]
syms h;                     % Vehicle altitude                          [m]
syms r;                     % Vehicle flight range                      [m]
syms t;                     % Vehicle time of flight                    [s]
syms g;                     % Acceleration                              [m/s^2]
syms Q;                     % Dynamic pressure at h                     [Pa]
syms W;                     % Vehicle weight                            [N]
syms C_D;                   % Drag Coefficient                          [adim]
syms A;                     % Vehicle reference area used in C_D        [m^2]
%syms beta;                    % Ballistic coefficient                     [N/m^2]
syms C_L;                   % Veicle lift coefficient                   [adim]

%% 2. PREVIOUS CALCULATIONS

% 2.1. Planet Input
Planet = getPlanet();
tic;

% 2.2. Define atmospheric variables
Rp = 0;             % Planet radius                         [m]
Mp = 0;             % Planet mass                           [kg]
g0 = 0;             % Acceleration at planet's surface      [m/s^2]
T0 = 0;             % US Standard Sea Level Temperature     [K]
P0 = 0;             % US Standard Sea Level Pressure        [Pa]
Mm = 0;             % Air molar mass                        [kg/mole]
H_layer = 0;        % Atmospheric layers altitude           [m]
lambda = 0;         % Atmospheric layers thermal gradients  [K/m]
gamma_gas = 0;      % Air specific heats relation           [adim]

if Planet == '1'        % Earth
    Rp = R_Earth;                   % Earth's radius                                [m]
    Mp = M_Earth;                   % Earth's mass                                  [kg]
    g0 = g0_Earth;                  % Earth's Surface Gravity Acceleration          [m/s^2]
    T0 = T0_Earth;                  % Earth's US Standard Sea Level Temperature     [K]
    P0 = P0_Earth;                  % Earth's US Standard Sea Level Pressure        [Pa]
    Mm = Mm_Earth;                  % Earth's air molar mass                        [kg/mole]
    H_layer = H_layer_Earth;        % Earth's atmospheric layers altitude           [m]
    lambda = lambda_layer_Earth;    % Earth's atmospheric layers thermal gradients  [K/m]
    gamma_gas = gamma_gas_Earth;    % Earth's air specific heats relation           [adim]
elseif Planet == '2'    % Mars
    Rp = R_Mars;                    % Mars's radius                                 [m]
    Mp = M_Mars;                    % Mars's mass                                   [kg]
    g0 = g0_Mars;                   % Mars's Surface Gravity Acceleration           [m/s^2]
    T0 = T0_Mars;                   % Mars's US Standard Sea Level Temperature      [K]
    P0 = P0_Mars;                   % Mars's US Standard Sea Level Pressure         [Pa]
    Mm = Mm_Mars;                   % Mars's air molar mass                         [kg/mole]
    H_layer = H_layer_Mars;         % Mars's atmospheric layers altitude            [m]
    lambda = lambda_layer_Mars;     % Mars's atmospheric layers thermal gradients   [K/m]
    gamma_gas = gamma_gas_Mars;     % Mars' air specific heats relation             [adim]
end

% 2.3. Compute base temperatures and pressures
[Tb, Pb] = getBaseTemperaturePressure(R, g0, T0, P0, Mm, H_layer, lambda);

plotAtmosphereV2(H_layer, lambda, Tb, Pb, R, g0, Mm);
if Planet == '1'
    
%% 3. SOLVER BALLISTIC REENTRY

% 3.1. Physical data
beta_ball = [100, 500, 1000, 5000];      % Balistic Coefficient [lbf/ft^2]
beta_ball = beta_ball*unitsratio('feet', 'meter')^2*4.4482216152605;  % Ballistic Coefficient [Pa]

% 3.2. Initial conditions
H0_ball = 250000*0.3048;         % Initial Altitude                      [m]
V0_ball = 22500*0.3048;          % Initial Velocity                      [m/s]     
gamma0_ball = 12;                % Initial Flight Path Angle             [deg]    

% 3.3. Numerical data
DeltaH = -10;

% 3.4. Declaration of arrays
H_ball = H0_ball:DeltaH:0;                                
V_ball = zeros(length(beta_ball), length(H_ball));
gamma_ball = zeros(length(beta_ball), length(H_ball));
t_ball = zeros(length(beta_ball), length(H_ball));
r_ball = zeros(length(beta_ball), length(H_ball));

a_ball = zeros(length(beta_ball), length(H_ball));
P_dynamic_ball = zeros(length(beta_ball), length(H_ball));
Ma_ball = zeros(length(beta_ball), length(H_ball));
Re_ball = zeros(length(beta_ball), length(H_ball));
P_stagnation_ball = zeros(length(beta_ball), length(H_ball));
h_stagnation_ball = zeros(length(beta_ball), length(H_ball));
HT_stagnation_ball = zeros(length(beta_ball), length(H_ball));
E_dynamic_ball = zeros(length(beta_ball), length(H_ball));

% 3.5. Computation of atmospheric properties at H_ball altitudes
[T_atmo, P_atmo, rho_atmo] = CalculateAtmosphere(Tb, Pb, H_layer, lambda, R, g0, Mm, H_ball);

% 3.6. Solve for every beta_ball
for i = 1:length(beta_ball)
    % Compute numerical solution using Runge-Kutta
    [V_ball(i,:), gamma_ball(i,:), t_ball(i,:), r_ball(i,:)] = ...
        numIntegration1_2(G, Mp, Rp, g0, R, Mm, Tb, Pb, H_layer, ...
        lambda, V0_ball, gamma0_ball, beta_ball(i), H_ball, DeltaH);
    % Compute other magnitudes
    a_ball(i,:) = computeAcceleration(V_ball(i,:), t_ball(i,:), g0_Earth);                  % Acceleration              [g0]
    P_dynamic_ball(i,:) = computeDynamicPressure(V_ball(i,:), rho_atmo);                    % Dynamic pressure          [Pa]
    Ma_ball(i,:) = computeMach(V_ball(i,:), T_atmo, Mm, R, gamma_gas);                      % Mach number               [adim]
    Re_ball(i,:) = computeReynolds(V_ball(i,:), T_atmo, rho_atmo);                          % Reynolds number           [adim]
    P_stagnation_ball(i,:) = ...
        computeStagnationPointPressure(V_ball(i,:), Ma_ball(i,:), gamma_gas, P_atmo);       % Stagnation pressure       [Pa]
    h_stagnation_ball(i,:) = computeStagnationPointEnthalpy(V_ball(i,:), T_atmo);           % Stagnation enthalpy       [J]
    HT_stagnation_ball(i,:) = computeStagnationPointHeatTransfer(V_ball(i,:), rho_atmo);    % Stagnation heat transfer  [W]
    E_dynamic_ball(i,:) = computeDynamicEnergy(V_ball(i,:), rho_atmo);                      % Dynamic energy            [J]
end

%% 4. SOLVER MERCURY CAPSULE

% 4.1. Physical data
W_merc = 2662.8;
S_merc = 30.27;
Length_merc = 6.2;
CD_merc = 1.60;
beta_merc = W_merc/(S_merc*CD_merc);
beta_merc = beta_merc*unitsratio('feet', 'meter')^2*4.4482216152605;

% 4.2. Initial conditions
H0_merc = 280000*0.3048;
V0_merc = 23000*0.3048;
gamma0_merc = 1.5;

% 4.3. Numerical data
DeltaH = -10;
H_merc = H0_merc:DeltaH:0;

% 4.4. Numerical integration
[V_merc, gamma_merc, t_merc, r_merc] = ...
    numIntegration1_2(G, Mp, Rp, g0, R, Mm, Tb, Pb, H_layer, lambda, ...
    V0_merc, gamma0_merc, beta_merc, H_merc, DeltaH);

% 4.5. Computation of other magnitudes
[T_atmo, P_atmo, rho_atmo] = CalculateAtmosphere(Tb, Pb, H_layer, lambda, R, g0, Mm, H_merc);
a_merc = computeAcceleration(V_merc, t_merc, g0_Earth);
Ma_merc = computeMach(V_merc, T_atmo, Mm, R, gamma_gas);

%% 5. SOLVER LIFTING ENTRY

% 5.1. Physical data
W_marv = 200000;                        % Weight                            [lbf]
S_marv = 2690;                          % Reference area                    [ft^2]
Length_marv = 107.5*0.3048;             % Reference lift                    [ft]
CL_marv = 0.84;                         % Lift coefficient                  [adim]
CD_marv = 0.84;                         % Drag coefficient                  [adim]
beta_marv = W_marv/(S_marv*CD_marv);	% Balistic coefficient              [lbf/ft^2]
beta_marv = beta_marv*unitsratio('feet', 'meter')^2*4.4482216152605; %      [Pa]

% 5.2. Initial conditions
H0_marv = 250000*0.3048;            % Initial Altitude              [m]
V0_marv = 23000*0.3048;             % Initial Velocity              [m/s]
gamma0_marv = [0.1 1.0 2.5];        % Initial Flight Path Angle     [deg]

% 5.3. Numerical data
t_max = 7200;   % Maximum time for numerical integration    [s]
Delta_t = 1;  % Time step                                   [s]

% 5.4. Declaration of arrays
t_marv = cell(length(gamma0_marv), 1);               % Time                      [s]
V_marv = cell(length(gamma0_marv), 1);               % Velocity                  [m/s]
gamma_marv = cell(length(gamma0_marv), 1);           % Flight path angle         [rad]
H_marv = cell(length(gamma0_marv), 1);               % Height                    [m]
r_marv = cell(length(gamma0_marv), 1);               % Range                     [m]
a_marv = cell(length(gamma0_marv), 1);               % Acceleration              [g0]
P_dynamic_marv = cell(length(gamma0_marv), 1);       % Dynamic pressure          [Pa]
Ma_marv = cell(length(gamma0_marv), 1);          % Mach number               [adim]
Re_marv = cell(length(gamma0_marv), 1);          % Reynolds number           [1/m]
P_stagnation_marv = cell(length(gamma0_marv), 1);    % Stagnation pressure       [Pa]
h_stagnation_marv = cell(length(gamma0_marv), 1);    % Stagnation enthalpy       [J/kg]
HT_stagnation_marv = cell(length(gamma0_marv), 1);   % Stagnation heat tranfer   [W/m^2]
E_dynamic_marv = cell(length(gamma0_marv), 1);       % Dynamic energy            [Pa�m/s]

% 5.5. Solve for every gamma0_marv
for i = 1:length(gamma0_marv)
    % Numerical integration
    [t_i, V_i, gamma_i, H_i, r_i] = ...
        numIntegration2_2(G, Mp, Rp, g0, R, Mm, Tb, Pb, H_layer, ...
        lambda, CL_marv, CD_marv, beta_marv, H0_marv, V0_marv, gamma0_marv(i), t_max, Delta_t);
    % Trim vectors to the length
    [t_i, V_i, gamma_i, H_i, r_i] = trimVectors(t_i, V_i, gamma_i, H_i, r_i);
    % Add to the respective matrices
    t_marv{i} = t_i;
    V_marv{i} = V_i;
    gamma_marv{i} = gamma_i;
    H_marv{i} = H_i;
    r_marv{i} = r_i;
    % Compute other magnitudes
    [T_atmo, P_atmo, rho_atmo] = ...
        CalculateAtmosphere(Tb, Pb, H_layer, lambda, R, g0, Mm, H_i);               % Atmospheric properties
    a_marv{i} = computeAcceleration(V_i, t_i, g0_Earth);                            % Acceleration              [g0]
    P_dynamic_marv{i} = computeDynamicPressure(V_i, rho_atmo);                      % Dynamic pressure          [Pa]
    Ma_marv{i} = computeMach(V_i, T_atmo, Mm, R, gamma_gas);                    % Mach number               [adim]
    Re_marv{i} = computeReynolds(V_i, T_atmo, rho_atmo).*Length_marv;           % Reynolds number           [adim]
    P_stagnation_marv{i} = ...
        computeStagnationPointPressure(V_i, Ma_marv{i}, gamma_gas, P_atmo);     % Stagnation pressure       [Pa]
    h_stagnation_marv{i} = computeStagnationPointEnthalpy(V_i, T_atmo);             % Stagnation enthalpy       [J]
    HT_stagnation_marv{i} = computeStagnationPointHeatTransfer(V_i, rho_atmo);      % Stagnation heat transfer  [W]
    E_dynamic_marv{i} = computeDynamicEnergy(V_i, rho_atmo);                        % Dynamic energy            [J]
end
elapsedTimeProgram = toc;
fprintf("%15s = %.5f %s", "Elapsed time", elapsedTimeProgram, "s");
%% 6. SOLUTION PLOTS

a = true;
fprintf("\n");
    while(a)
        disp('Insert the solution which you want to visualise');
        disp('You can choose between ballistic, mercury or lifting entry');
        prompt = '[Answer (B(Ballistic)/M(Mercury)/L(Lifting)/A(All of them)]:';
        solution = input(prompt,'s');
        
        % Get plots from user input
        if ((solution == "B") || (solution == "Ballistic") || (solution == "b"))
            a = false;
            
            % 3.7. Plot ballistic solution
            plotSolutionBallistic(beta_ball, H_ball, V_ball, gamma_ball, t_ball, r_ball, a_ball, ...
              P_dynamic_ball, Ma_ball, Re_ball, P_stagnation_ball, h_stagnation_ball, ...
              HT_stagnation_ball, E_dynamic_ball, 1);
        
        elseif ((solution == "M") || (solution == "Mercury") || (solution == "m"))
            a = false;
            
            % 4.6. Plot mercury solution
            plotSolutionMercury(H_merc, V_merc, t_merc, r_merc, a_merc, Ma_merc, 1);
            
        elseif ((solution == "L") || (solution == "Lifting") || (solution == "l"))
            a = false;
            
            % 5.6. Plot solution
            plotSolutionMARV(gamma0_marv, H_marv, V_marv, gamma_marv, t_marv, r_marv, a_marv, ...
               P_dynamic_marv, Ma_marv, Re_marv, P_stagnation_marv, ...
               h_stagnation_marv, HT_stagnation_marv, E_dynamic_marv, 1);
        elseif ((solution == "A") || (solution == "All of them") || (solution == "a"))
            a = false;
            
            % 3.7. Plot ballistic solution
            plotSolutionBallistic(beta_ball, H_ball, V_ball, gamma_ball, t_ball, r_ball, a_ball, ...
              P_dynamic_ball, Ma_ball, Re_ball, P_stagnation_ball, h_stagnation_ball, ...
              HT_stagnation_ball, E_dynamic_ball, 1);
            % 4.6. Plot mercury solution
            plotSolutionMercury(H_merc, V_merc, t_merc, r_merc, a_merc, Ma_merc, 1);
            % 5.6. Plot solution
            plotSolutionMARV(gamma0_marv, H_marv, V_marv, gamma_marv, t_marv, r_marv, a_marv, ...
               P_dynamic_marv, Ma_marv, Re_marv, P_stagnation_marv, ...
               h_stagnation_marv, HT_stagnation_marv, E_dynamic_marv, 1);
        else
            disp('    ')
            disp('Invalid answer. Please re-enter the answer.')
            disp('    ')
        end
    end
elapsedTimePlots = toc - elapsedTimeProgram;
fprintf("%15s = %.5f %s", "Plotting Elapsed time", elapsedTimePlots, "s");

end






