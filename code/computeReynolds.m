function Re = computeReynolds(V, T, rho)
%--------------------------------------------------------------------------
% Inputs:
%   - V         Temperature [K]
%   - T         Temperature [K]
%   - rho       Temperature [K]
%--------------------------------------------------------------------------
% Outputs:
%   - Re        Reynolds number [1/m]
% Notes:
%   - There is no reference longitude, so Re number is not 
%--------------------------------------------------------------------------

Re = zeros(1, length(V));
for i = 1:length(V)
    Re(i) = rho(i)*V(i)/computeDynamicViscosity(T(i));
end

end