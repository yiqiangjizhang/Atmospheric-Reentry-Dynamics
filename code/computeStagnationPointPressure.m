function P_stagnation = computeStagnationPointPressure(V, Ma, gamma_gas, P_atmo)
%--------------------------------------------------------------------------
% Inputs:
%   - V                 Vector of velocities                        [m/s]
%   - Ma                Vector of mach numbers                      [adim]
%   - gamma_gas         Relation of specific heats of air           [adim] 
%   - P_atmo            Atmospheric pressure at solution altitudes  [Pa]   
%--------------------------------------------------------------------------
% Outputs:
%   - P_stagnation      Stagnation point pressure  [Pa]
%--------------------------------------------------------------------------

P_stagnation = zeros(1, length(V));
for i = 1:length(V)
    if Ma(i) < 1
        P_stagnation(i) = P_atmo(i)*(1 + ((gamma_gas-1)/2)*Ma(i)^2)^(gamma_gas/(gamma_gas-1));
    elseif Ma(i) >= 1
        P_shock = P_atmo(i)*(1 + (2*gamma_gas)/(gamma_gas+1)*(Ma(i)^2-1));
        P_stagnation(i) = P_shock*(1 + ((gamma_gas-1)/2)*(1 + ((gamma_gas-1)/2)*Ma(i)^2)/(gamma_gas*Ma(i)^2-(gamma_gas-1)/2))^(gamma_gas/(gamma_gas-1));
    end
end

end

