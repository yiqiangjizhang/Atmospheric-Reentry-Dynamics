function P = getPressureV2(Tb, Pb, H_layer, lambda, R, g0, Mm, H)
%--------------------------------------------------------------------------
% Inputs:
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - Pb        Vector of base pressures, 1 x length(H_layer)       [Pa]
%   - H_layer   Atmospheric layers altitude                         [m]
%   - lambda    Atmospheric layers thermal gradients                [K/m]
%   - R         Universal gas constant                              [J/(kg K)]
%   - g0        Acceleration at planet's surface                    [m/s^2]
%   - Mm        Planet's air molar mass                             [kg/mol]
%   - H         Altitude                                            [m]
%--------------------------------------------------------------------------
% Outputs:
%   - P         Pressure at altitude H
%--------------------------------------------------------------------------


P = 0;
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
    if lambda(layer) == 0
        P = Pb(layer)*exp(-g0*Mm*(H-H_layer(layer))/(R*Tb(layer)));
    else
        P = Pb(layer)*((Tb(layer)/(Tb(layer) + lambda(layer)*(H-H_layer(layer))))^(g0*Mm/(R*lambda(layer))));
    end
end

end

