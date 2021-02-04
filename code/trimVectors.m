function [t, V, gamma, H, r] = trimVectors(t, V, gamma, H, r)
%--------------------------------------------------------------------------
% Inputs:
%       - t             Time vector                 [s]
%       - V             Velocities vector           [m/s]
%       - gamma         Flight path angle vector    [rad]
%       - H             Height vector               [m]
%       - r             Range vector                [m]
% Note: when doin the MARV numerical integration, the final condition is
% H=0. This condition is achieved before time reaches its maximum, so
% numerical integration has to stop. Therefore, approximately half of the
% solution vector is empty and must be erased.
%--------------------------------------------------------------------------
% Outputs:
%       - t             Time vector                 [s]
%       - V             Velocities vector           [m/s]
%       - gamma         Flight path angle vector    [rad]
%       - H             Height vector               [m]
%       - r             Range vector                [m]
%--------------------------------------------------------------------------

% Find which is the last index with physical sense (H(i) > 0)
i = 1;      % Index
found = 0;  % Condition
while (i <= length(H)) && found == 0
    if H(i) < 0
        found = 1;
    else
        i = i + 1;
    end
end

% Trim vectors
t(i:end) = [];
V(i:end) = [];
gamma(i:end) = [];
H(i:end) = [];
r(i:end) = [];

end

