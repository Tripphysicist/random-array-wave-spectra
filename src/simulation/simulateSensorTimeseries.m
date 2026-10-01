function [eta, Phi_flat] = simulateSensorTimeseries(Spec, sensorX, sensorY, t, Phi_flat)
% simulateSensorTimeseries Simulates wave elevation at specific locations.
% Optionally accepts a fixed phase array to sample the same wave field
% with different sensor geometries.

g = 9.81; % gravity
w = Spec.w(:);
theta = Spec.theta(:);

S = Spec.S; 
if size(S, 1) == numel(theta) && size(S, 2) == numel(w)
    S = S'; 
end

dw = gradient(w); 
dtheta = gradient(theta);

[W, Theta] = ndgrid(w, theta);
[dW, dTheta] = ndgrid(dw, dtheta);

K = (W.^2) / g; 
Kx = K .* cos(Theta);
Ky = K .* sin(Theta);

A = sqrt(2 * S .* dW .* dTheta);

A_flat = A(:);
Kx_flat = Kx(:);
Ky_flat = Ky(:);
W_flat = W(:);

% Filter near-zero energy components to increase speed
idx = A_flat > 1e-6 * max(A_flat);
A_flat = A_flat(idx);
Kx_flat = Kx_flat(idx);
Ky_flat = Ky_flat(idx);
W_flat = W_flat(idx);
numComponents = numel(A_flat);

% If phases are not provided, generate new ones
if nargin < 5 || isempty(Phi_flat)
    Phi_flat = 2 * pi * rand(numComponents, 1);
end

numSensors = numel(sensorX);
numTimeSteps = numel(t);
eta = zeros(numTimeSteps, numSensors);
t_col = t(:); 

% PRE-COMPUTE TIME MATRIX OUTSIDE SENSOR LOOP
% W_flat' is [1 x Nc], t_col is [Nt x 1]
% Wt_mat becomes an [Nt x Nc] matrix
Wt_mat = t_col * W_flat'; 
cos_Wt = cos(Wt_mat);
sin_Wt = sin(Wt_mat);

for s = 1:numSensors
    % Calculate spatial phase for this sensor [Nc x 1]
    spatialPhase = Kx_flat * sensorX(s) + Ky_flat * sensorY(s);
    totalPhase = spatialPhase + Phi_flat;
    
    % Amplitude and Phase vectors [Nc x 1]
    A_cos_phase = A_flat .* cos(totalPhase); 
    A_sin_phase = A_flat .* sin(totalPhase); 
    
    % Matrix multiplication intrinsically sums across the Nc dimension!
    % [Nt x Nc] * [Nc x 1] = [Nt x 1] time series vector
    eta(:, s) = cos_Wt * A_cos_phase + sin_Wt * A_sin_phase;
end
end