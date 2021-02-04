function [Tb, Pb] = getBaseTemperaturePressure(R, g0, T0, P0, Mm, H_layer, lambda)
%--------------------------------------------------------------------------
% Inputs:
%   - R         Universal Constant for Ideal Gases          [J/mole*K]
%   - g0        Acceleration at planet's surface            [m/s^2]
%   - T0        Standard Temperature at planet's surface    [K]
%   - P0        Standard Pressure at planet's surface       [Pa]
%   - Mm        Atmospheric gas molecular mass              [kg*mole^-1]    
%   - H_layer   Atmospheric layers altitude                 [m]
%   - lambda    Atmospheric layers thermal gradients        [K/m]
%--------------------------------------------------------------------------
% Outputs:
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - Pb        Vector of base pressures, 1 x length(H_layer)       [Pa]
% Base (magnitude) refers to the value of (magnitude( at the beginning of
% an atmospheric layer
%--------------------------------------------------------------------------

% Declare vectors of Base Temperature and Pressure
Tb = zeros(1, length(H_layer)); % Base Temperature  [K]
Pb = zeros(1, length(H_layer)); % Base Pressure     [Pa]

% Compute Base Temperatures and Pressures for each layer
Tb(1) = T0;
Pb(1) = P0;
for i = 2:length(H_layer)
    % Compute Base Temperature at layer i
    Tb(i) = Tb(i-1) + lambda(i-1)*(H_layer(i)-H_layer(i-1));
    % Compute Base Pressure at layer i
    if lambda(i-1) == 0 % Isothermal layer
        Pb(i) = Pb(i-1)*exp(-g0*Mm*(H_layer(i)-H_layer(i-1))/(R*Tb(i-1)));
    else                % Non-isothermal layer
        Pb(i) = Pb(i-1)*((Tb(i-1)/(Tb(i-1)+lambda(i-1)*(H_layer(i)-H_layer(i-1))))^(g0*Mm/(R*lambda(i-1))));
    end
end

% fprintf("%8s%12s%12s%12s\n", "Layer", "Height(m)", "T(K)", "P(Pa)");
% for i = 1:length(H_layer)
%     fprintf("%8d", i);
%     fprintf("%12d", H_layer(i));
%     fprintf("%12.2f", Tb(i));
%     fprintf("%12.2f\n", Pb(i));
% end
% fprintf("\n\n");

end