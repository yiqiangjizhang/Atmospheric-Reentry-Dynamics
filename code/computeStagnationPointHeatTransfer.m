function HT_stagnation = computeStagnationPointHeatTransfer(V_sol, rho)
%--------------------------------------------------------------------------
% Inputs:
%   - V_sol             Vector of velocities            [m/s]
%   - rho               Vector of densities             [kg/m^3]
%--------------------------------------------------------------------------
% Outputs:
%   - HT_stagnation     Stagnation point heat transfer  [W/m^2]
%--------------------------------------------------------------------------

%k_Earth_NASA = 1.7415*10^-4;
k_Earth = 1.83 * 10^-4;
%k_Mars_NASA = 1.9027*10^-4;
k_Mars = 1.89 * 10^-4;
R_n = 1*0.3048; %bluntness sphere

HT_stagnation = -k_Earth*sqrt(rho/R_n).*V_sol.^3;

end