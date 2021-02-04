function h_stagnation = computeStagnationPointEnthalpy(V_sol, T_sol)
%--------------------------------------------------------------------------
% Inputs:
%   - V_sol             Vector of velocities            [m/s]
%   - T_sol             Vector of temperatues           [K]
%--------------------------------------------------------------------------
% Outputs:
%   - h_stagnation      Stagnation point enthalpy       [J/kg]
%--------------------------------------------------------------------------

C_p = 1034.9 - 2.849e-1*T_sol +7.817e-4*T_sol.^2 -4.971e-7*T_sol.^3 + 1.007e-10*T_sol.^4;
h_stagnation =  C_p.*T_sol + 1/2*V_sol.^2;


end