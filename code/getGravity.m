function g = getGravity(H, G, R, M)
%--------------------------------------------------------------------------
% Inputs:
%   - H        Vector of altitudes, 1 x length(H_layer)            [m]
%   - G        Gravitational Constant                              [m^3/(kg*s^2)]
%   - R        Universal Constant for Ideal Gases                  [J/mole*K]
%   - M        Molecular Mass                                      [kg*mole^-1]                         [m]
%--------------------------------------------------------------------------
% Outputs:
%   - g        Vector of gravities, 1 x length(H_layer)            [m/s^2]
% Base (magnitude) refers to the value of (magnitude( at the beginning of
% an atmospheric layer
%--------------------------------------------------------------------------

g = G*M/(R+H)^2;    % [m/s^2]

end