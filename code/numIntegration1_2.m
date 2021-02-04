function [V_sol, gamma_sol, t_sol, r_sol] = numIntegration1_2(G, Mp, Rp, g0, R, Mm, Tb, Pb, H_layer, lambda, V0, gamma0, beta, H_sol, DeltaH)
%--------------------------------------------------------------------------
% Inputs:
%   - Constants:
%       - G         Gravitational Constant                              [m^3/(kg*s^2)]
%       - Mp        Planet mass                                         [kg]
%       - Rp        Planet radius                                       [m]
%       - g0        Acceleration at planet's surface                    [m/s^2]
%       - R         Universal Constant for Ideal Gases                  [J/mole*K]
%       - Mm        Air molar mass                                      [kg/mole]
%       - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%       - Pb        Vector of base pressures, 1 x length(H_layer)       [K]
%       - H_layer   Atmospheric layers altitude                         [m]
%       - lambda    Atmospheric layers thermal gradients                [K/m]
%   - Physical data:
%       - beta      Balistic coefficient                                [Pa]
%   - Initial conditions:
%       - H0        Initial altitude                                    [m]
%       - V0        Initial velocity                                    [m/s]
%       - gamma     Initial Flight Path Angle                           [deg]
%   - Numerical data:
%       - DeltaH    Altitude step                                       [m]
%--------------------------------------------------------------------------
% Outputs: solution to the ODE system
%   - H_sol         Altitudes vector                                    [m]
%   - V_sol         Velocities vector                                   [m/s]
%   - gamma_sol     Flight path angles vector                           [rad]
%   - t_sol         Time vector                                         [s]
%   - r_sol         Ranges vector                                       [m]
%--------------------------------------------------------------------------

% Declaration of solution vectors
V_sol = zeros(1, length(H_sol));
gamma_sol = zeros(1, length(H_sol));
t_sol = zeros(1, length(H_sol));
r_sol = zeros(1, length(H_sol));

% Initial conditions
V_sol(1) = V0;
gamma_sol(1) = deg2rad(gamma0);

% Runge-Kutta
for i = 1:length(H_sol)-1
    % Definition of functions
    F1 = @(H, V, gamma, t, r) (G*Mp/(Rp+H)^2)*(getDensityV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H)*V^2/(2*beta)  - sin(gamma))/(V*sin(gamma));
    F2 = @(H, V, gamma, t, r) ((V^2/(Rp+H) - G*Mp/(Rp+H)^2)/V^2)*(cos(gamma)/sin(gamma));
    F3 = @(H, V, gamma, t, r) -1/(V*sin(gamma));
    F4 = @(H, V, gamma, t, r) -(Rp/(Rp+H))*cos(gamma)/sin(gamma);
    % Computation of coefficients sub 1
    i1 = F1(H_sol(i), V_sol(i), gamma_sol(i), t_sol(i), r_sol(i));
    j1 = F2(H_sol(i), V_sol(i), gamma_sol(i), t_sol(i), r_sol(i));
    k1 = F3(H_sol(i), V_sol(i), gamma_sol(i), t_sol(i), r_sol(i));
    l1 = F4(H_sol(i), V_sol(i), gamma_sol(i), t_sol(i), r_sol(i));
    % Computation of coefficients sub 2
    i2 = F1(H_sol(i)+DeltaH/2, V_sol(i)+i1*DeltaH/2, gamma_sol(i)+j1*DeltaH/2, t_sol(i)+k1*DeltaH/2, r_sol(i)+l1*DeltaH/2);
    j2 = F2(H_sol(i)+DeltaH/2, V_sol(i)+i1*DeltaH/2, gamma_sol(i)+j1*DeltaH/2, t_sol(i)+k1*DeltaH/2, r_sol(i)+l1*DeltaH/2);
    k2 = F3(H_sol(i)+DeltaH/2, V_sol(i)+i1*DeltaH/2, gamma_sol(i)+j1*DeltaH/2, t_sol(i)+k1*DeltaH/2, r_sol(i)+l1*DeltaH/2);
    l2 = F4(H_sol(i)+DeltaH/2, V_sol(i)+i1*DeltaH/2, gamma_sol(i)+j1*DeltaH/2, t_sol(i)+k1*DeltaH/2, r_sol(i)+l1*DeltaH/2);
    % Computation of coefficients sub 3
    i3 = F1(H_sol(i)+DeltaH/2, V_sol(i)+i2*DeltaH/2, gamma_sol(i)+j2*DeltaH/2, t_sol(i)+k2*DeltaH/2, r_sol(i)+l2*DeltaH/2);
    j3 = F2(H_sol(i)+DeltaH/2, V_sol(i)+i2*DeltaH/2, gamma_sol(i)+j2*DeltaH/2, t_sol(i)+k2*DeltaH/2, r_sol(i)+l2*DeltaH/2);
    k3 = F3(H_sol(i)+DeltaH/2, V_sol(i)+i2*DeltaH/2, gamma_sol(i)+j2*DeltaH/2, t_sol(i)+k2*DeltaH/2, r_sol(i)+l2*DeltaH/2);
    l3 = F4(H_sol(i)+DeltaH/2, V_sol(i)+i2*DeltaH/2, gamma_sol(i)+j2*DeltaH/2, t_sol(i)+k2*DeltaH/2, r_sol(i)+l2*DeltaH/2);
    % Computation of coefficients sub 4
    i4 = F1(H_sol(i)+DeltaH, V_sol(i)+i3*DeltaH, gamma_sol(i)+j3*DeltaH, t_sol(i)+k3*DeltaH, r_sol(i)+l3*DeltaH);
    j4 = F2(H_sol(i)+DeltaH, V_sol(i)+i3*DeltaH, gamma_sol(i)+j3*DeltaH, t_sol(i)+k3*DeltaH, r_sol(i)+l3*DeltaH);
    k4 = F3(H_sol(i)+DeltaH, V_sol(i)+i3*DeltaH, gamma_sol(i)+j3*DeltaH, t_sol(i)+k3*DeltaH, r_sol(i)+l3*DeltaH);
    l4 = F4(H_sol(i)+DeltaH, V_sol(i)+i3*DeltaH, gamma_sol(i)+j3*DeltaH, t_sol(i)+k3*DeltaH, r_sol(i)+l3*DeltaH);
    % Compute next step
    V_sol(i+1) = V_sol(i) + (DeltaH/6)*(i1 + 2*i2 + 2*i3 + i4);
    gamma_sol(i+1) = gamma_sol(i) + (DeltaH/6)*(j1 + 2*j2 + 2*j3 + j4);
    t_sol(i+1) = t_sol(i) + (DeltaH/6)*(k1 + 2*k2 + 2*k3 + k4);
    r_sol(i+1) = r_sol(i) + (DeltaH/6)*(l1 + 2*l2 + 2*l3 + l4);
end


end

