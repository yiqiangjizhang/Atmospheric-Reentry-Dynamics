function E_dynamic = computeDynamicEnergy(V, rho)
%--------------------------------------------------------------------------
% Inputs:
%   - V             Velocity vector         [m/s]
%   - rho           Densities vector        [s]
%--------------------------------------------------------------------------
% Outputs: 
%   - E_dynamic     Dynamic energy vector   [Pa·m/s]
%--------------------------------------------------------------------------

E_dynamic = 1/2 * rho .* V.^3; 

end