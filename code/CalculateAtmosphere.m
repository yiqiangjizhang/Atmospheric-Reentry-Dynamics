function [T_sol, P_sol, rho_sol] = CalculateAtmosphere(Tb, Pb, H_layer, lambda, R, g0, Mm, H_sol)
%--------------------------------------------------------------------------
% Inputs:
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - H_layer   Atmospheric layers altitude                         [m]
%   - lambda    Atmospheric layers thermal gradients                [K/m]
%   - R         Universal gas constant                              [J/(kg K)]
%   - g0        Acceleration at planet's surface                    [m/s^2]
%   - Mm        Planet's air molar mass                             [kg/mol]
%   - H_sol     Vector of altitudes to compute atmosphere           [m]
%--------------------------------------------------------------------------
% Outputs:
%   - T_sol     Temperatures at altitudes defined by H_sol
%   - P_sol     Pressures at altitudes defined by H_sol
%   - rho_sol   Densities at altitudes defined by H_sol
%--------------------------------------------------------------------------


T_sol = zeros(1, length(H_sol));
P_sol = zeros(1, length(H_sol));
rho_sol = zeros(1, length(H_sol));

for i = 1:length(H_sol)
    T_sol(i) = getTemperatureV2(Tb, H_layer, lambda, H_sol(i));
    P_sol(i) = getPressureV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H_sol(i));
    rho_sol(i) = getDensityV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H_sol(i));
end

end

