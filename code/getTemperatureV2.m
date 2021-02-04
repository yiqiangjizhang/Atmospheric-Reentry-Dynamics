function T = getTemperatureV2(Tb, H_layer, lambda, H)
%--------------------------------------------------------------------------
% Inputs:
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - H_layer   Atmospheric layers altitude                         [m]
%   - lambda    Atmospheric layers thermal gradients                [K/m]
%   - H         Altitude at which 
%--------------------------------------------------------------------------
% Outputs:
%   - Pb        Vector of base pressures, 1 x length(H_layer)       [Pa]
% Base (magnitude) refers to the value of (magnitude( at the beginning of
% an atmospheric layer
%--------------------------------------------------------------------------

T = 0;
found = 0;
layer = 1;

while (layer <= length(H_layer)-1) && (found == 0)
    if (H_layer(layer) <= H) && (H < H_layer(layer+1))
        found = 1;
    else
        layer = layer + 1;
    end
end

if found == 1
    T = Tb(layer) + lambda(layer)*(H - H_layer(layer)); 
end

end

