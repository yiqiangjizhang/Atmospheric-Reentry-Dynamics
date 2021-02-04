function Ma = computeMach(V, T, Mm, R, gamma_gas)
%--------------------------------------------------------------------------
% Inputs:
%   - V         Temperature [K]
%   - T         Temperature [K]
%   - rho       Temperature [K]
%--------------------------------------------------------------------------
% Outputs:
%   - Ma        Mach number [adim]
%--------------------------------------------------------------------------

Ma = V./(sqrt(gamma_gas*(R/Mm)*T));

end