function [t_sol, V_sol, gamma_sol, H_sol, r_sol] = ...
    numIntegration2_2(G, Mp, Rp, g0, R, Mm, Tb, Pb, H_layer, lambda, ...
    CL, CD, beta, H0, V0, gamma0, t_max, Delta_t)
%--------------------------------------------------------------------------
% Inputs:
%   - Constants:
%       - G             Gravitational Constant                              [m^3/(kg*s^2)]
%       - Mp            Planet mass                                         [kg]
%       - Rp            Planet radius                                       [m]
%       - g0            Acceleration at planet's surface                    [m/s^2]
%       - R             Universal Constant for Ideal Gases                  [J/mole*K]
%       - Mm            Air molar mass                                      [kg/mole]
%       - Tb            Vector of base temperatures, 1 x length(H_layer)    [K]
%       - Pb            Vector of base pressures, 1 x length(H_layer)       [K]
%       - H_layer       Atmospheric layers altitude                         [m]
%       - lambda        Atmospheric layers thermal gradients                [K/m]
%   - Physical data:
%       - CL            Lift coefficient            [adim]
%       - CD            Drag coefficient            [adim]
%       - beta          Balistic coefficient        [Pa]
%   - Initial conditions:
%       - H0            Initial altitude            [m]
%       - V0            Initial velocity            [m/s]
%       - gamma0        Initial Flight Path Angle   [deg]
%   - Numerical data:
%       - t_max         Max time                    [s]
%       - Delta_t       Time step                   [s]
%--------------------------------------------------------------------------
% Outputs: solution to the ODE system
%   - t_sol             Time vector                 [s]
%   - H_sol             Altitudes vector            [m]
%   - V_sol             Velocities vector           [m/s]
%   - gamma_sol         Flight path angles vector	[rad]
%   - r_sol             Range vectors               [m]
%--------------------------------------------------------------------------


% Declaration of solution vectors
t_sol = 0:Delta_t:t_max;
H_sol = zeros(1, length(t_sol));
V_sol = zeros(1, length(t_sol));
gamma_sol = zeros(1, length(t_sol));
r_sol = zeros(1, length(t_sol));

% Initial conditions
H_sol(1) = H0;
V_sol(1) = V0;
gamma_sol(1) = deg2rad(gamma0);

% Runge-Kutta
i = 1;
while (i <= length(t_sol)-1) && (H_sol(i) >= 0)
    % Definition of functions
    F1 = @(t, V, gamma, H, r) ((G*Mp)/(Rp+H)^2)*(sin(gamma)-((1/2)*getDensityV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H)*V^2)/beta);
    F2 = @(t, V, gamma, H, r) ((-(1/2)*getDensityV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H)*V^2)*((G*Mp)/(Rp+H)^2)*(1/beta)*(CL/CD) + (((G*Mp)/(Rp+H)^2)-(V^2/(Rp+H)))*cos(gamma))/V;
    F3 = @(t, V, gamma, H, r) -V*sin(gamma);
    F4 = @(t, V, gamma, H, r) (Rp/(Rp+H))*V*cos(gamma);
    % Computation of coefficients sub 1
    i1 = F1(t_sol(i), V_sol(i), gamma_sol(i), H_sol(i), r_sol(i));
    j1 = F2(t_sol(i), V_sol(i), gamma_sol(i), H_sol(i), r_sol(i));
    k1 = F3(t_sol(i), V_sol(i), gamma_sol(i), H_sol(i), r_sol(i));
    l1 = F4(t_sol(i), V_sol(i), gamma_sol(i), H_sol(i), r_sol(i));    
    % Computation of coefficients sub 2
    i2 = F1(t_sol(i)+Delta_t/2, V_sol(i)+i1*Delta_t/2, gamma_sol(i)+j1*Delta_t/2, H_sol(i)+k1*Delta_t/2, r_sol(i)+l1*Delta_t/2);
    j2 = F2(t_sol(i)+Delta_t/2, V_sol(i)+i1*Delta_t/2, gamma_sol(i)+j1*Delta_t/2, H_sol(i)+k1*Delta_t/2, r_sol(i)+l1*Delta_t/2);
    k2 = F3(t_sol(i)+Delta_t/2, V_sol(i)+i1*Delta_t/2, gamma_sol(i)+j1*Delta_t/2, H_sol(i)+k1*Delta_t/2, r_sol(i)+l1*Delta_t/2);
    l2 = F4(t_sol(i)+Delta_t/2, V_sol(i)+i1*Delta_t/2, gamma_sol(i)+j1*Delta_t/2, H_sol(i)+k1*Delta_t/2, r_sol(i)+l1*Delta_t/2);
    % Computation of coefficients sub 3
    i3 = F1(t_sol(i)+Delta_t/2, V_sol(i)+i2*Delta_t/2, gamma_sol(i)+j2*Delta_t/2, H_sol(i)+k2*Delta_t/2, r_sol(i)+l2*Delta_t/2);
    j3 = F2(t_sol(i)+Delta_t/2, V_sol(i)+i2*Delta_t/2, gamma_sol(i)+j2*Delta_t/2, H_sol(i)+k2*Delta_t/2, r_sol(i)+l2*Delta_t/2);
    k3 = F3(t_sol(i)+Delta_t/2, V_sol(i)+i2*Delta_t/2, gamma_sol(i)+j2*Delta_t/2, H_sol(i)+k2*Delta_t/2, r_sol(i)+l2*Delta_t/2);
    l3 = F4(t_sol(i)+Delta_t/2, V_sol(i)+i2*Delta_t/2, gamma_sol(i)+j2*Delta_t/2, H_sol(i)+k2*Delta_t/2, r_sol(i)+l2*Delta_t/2);
    % Computation of coefficients sub 4
    i4 = F1(t_sol(i)+Delta_t, V_sol(i)+i3*Delta_t, gamma_sol(i)+j3*Delta_t, H_sol(i)+k3*Delta_t, r_sol(i)+l3*Delta_t);
    j4 = F2(t_sol(i)+Delta_t, V_sol(i)+i3*Delta_t, gamma_sol(i)+j3*Delta_t, H_sol(i)+k3*Delta_t, r_sol(i)+l3*Delta_t);
    k4 = F3(t_sol(i)+Delta_t, V_sol(i)+i3*Delta_t, gamma_sol(i)+j3*Delta_t, H_sol(i)+k3*Delta_t, r_sol(i)+l3*Delta_t);
    l4 = F4(t_sol(i)+Delta_t, V_sol(i)+i3*Delta_t, gamma_sol(i)+j3*Delta_t, H_sol(i)+k3*Delta_t, r_sol(i)+l3*Delta_t);
    % Compute next step
    V_sol(i+1) = V_sol(i) + (Delta_t/6)*(i1 + 2*i2 + 2*i3 + i4);
    gamma_sol(i+1) = gamma_sol(i) + (Delta_t/6)*(j1 + 2*j2 + 2*j3 + j4);
    H_sol(i+1) = H_sol(i) + (Delta_t/6)*(k1 + 2*k2 + 2*k3 + k4);
    r_sol(i+1) = r_sol(i) + (Delta_t/6)*(l1 + 2*l2 + 2*l3 + l4);
    i = i + 1;
end

end

