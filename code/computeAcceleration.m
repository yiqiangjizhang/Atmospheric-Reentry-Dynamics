function a = computeAcceleration(V, t, g0_Earth)
%--------------------------------------------------------------------------
% Inputs:
%   - V         Velocity vector                     [m/s]
%   - t         Time vector                         [s]
%   - g0_Earth  Acceleration at Earth's surface     [m/s^2]
%--------------------------------------------------------------------------
% Outputs: the function computes the acceleration at every point of the
% trajectory using a 2nd order scheme
%   - a         Acceleration vector                 [g0]
%--------------------------------------------------------------------------

% Computation of acceleration
a = zeros(1, length(V));
for i = 2:length(V)-1
    a(i) = (V(i+1)-V(i-1))/(t(i+1)-t(i-1));
end

% Normalization using g0_Earth
a = -a;
a(1) = a(2);
a(end-1) = a(end);
a = a/g0_Earth;

end

