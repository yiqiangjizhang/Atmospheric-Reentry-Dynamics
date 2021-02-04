function P_dynamic = computeDynamicPressure(V, rho)
%--------------------------------------------------------------------------
% Inputs:
%   - V             Vector of velocities    [m/s]
%   - rho           Vector of densities     [K]
%--------------------------------------------------------------------------
% Outputs:
%   - P_dynamic     Dynamic pressure        [Pa]
%--------------------------------------------------------------------------

P_dynamic = (1/2)*rho.*(V.^2);


end
