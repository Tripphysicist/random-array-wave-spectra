% --- Sensor positions (2D coordinates) ---


 xp = x(:,end,I(1));
 yp = y(:,end,I(1));

% Example: random array in a 10m x 10m box
N = length(xp);
X = [xp(:), yp(:)];  % Ensure column vectors and combine

% --- Compute co-array (difference vectors between sensor pairs) ---
coarray = [];
for i = 1:N
    for j = 1:N
        r = X(i,:) - X(j,:);
        coarray(end+1, :) = r;
    end
end

% --- Convert to wavenumber space (kx, ky) ---
% Normalize each coordinate by dividing by wavelength (or multiply by wavenumber)
% We'll scale by kp to examine resolution at peak frequency
fp = 1/6;  % Hz (peak frequency)
g = 9.81;
kp = (2*pi*fp)^2 / g;

% Scale to dimensionless wavenumber coordinates
kx = kp * coarray(:,1);
ky = kp * coarray(:,2);

% --- Plot ---
figure;
scatter(kx, ky, 10, 'filled');
xlabel('k_x / k_p');
ylabel('k_y / k_p');
axis equal;
title('Normalized Co-array Coverage');
grid on;
hold on

% --- Overlay unit circle (corresponds to kp) ---
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), '--k', 'LineWidth', 1.2);
text(cos(pi/4), sin(pi/4), 'k = k_p', 'FontSize', 8);

%% array response

% Parameters
N = length(x);
k_max = 2;               % Range of wavenumber space to explore (in units of k_p)
res = 200;               % Resolution of the grid

% Grid in wavenumber space (kx, ky)
kx = linspace(-k_max, k_max, res);
ky = linspace(-k_max, k_max, res);
[KX, KY] = meshgrid(kx, ky);

% Compute array response
ARF = zeros(size(KX));
for n = 1:N
    ARF = ARF + exp(1i * (KX * x(n) + KY * y(n)));
end
ARF = abs(ARF).^2;  % Power response

% Normalize and plot
ARF = ARF / max(ARF(:));  % Normalize
figure;
imagesc(kx, ky, ARF);
axis xy equal;
xlabel('$k_x / k_p$', 'interpreter', 'latex');
ylabel('$k_y / k_p$', 'interpreter', 'latex');
title('Array Response Function');
colorbar;
