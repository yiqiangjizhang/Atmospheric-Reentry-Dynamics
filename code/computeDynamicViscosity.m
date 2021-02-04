function mu = computeDynamicViscosity(T)
%--------------------------------------------------------------------------
% Inputs:
%   - T     Temperature [K]
%--------------------------------------------------------------------------
% Output: 
%   - mu    Dynamic viscosity at temperature T  [Pa·s]
%--------------------------------------------------------------------------

B = 1.458e-6;               % Constant              [kg/(s m K^(1/2))]
S = 110.4;                  % Sutherland's constant [K]
mu = (B*T^(3/2))/(T + S);   % Dynamic viscosity     [Pa s]

end