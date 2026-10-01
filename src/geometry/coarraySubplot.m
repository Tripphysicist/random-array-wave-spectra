function coarraySubplot(xp, yp, Tp, color)
    % xp, yp: [N x 1] sensor coordinates
    % Tp: peak period
    % color: RGB vector for the specific configuration
    
    % --- Physics ---
    fp = 1/Tp;
    g = 9.81;
    kp = (2*pi*fp)^2 / g; 

    % --- Compute Co-array (Vectorized) ---
    N = length(xp);
    [idx1, idx2] = meshgrid(1:N, 1:N);
    % Scaling: kp [rad/m] * spatial_lag [m] = dimensionless phase space
    kx_lag = kp * (xp(idx1(:)) - xp(idx2(:)));
    ky_lag = kp * (yp(idx1(:)) - yp(idx2(:)));

    % --- Plot into existing subplot ---
    hold on;
    plot(kx_lag, ky_lag, '.', 'MarkerSize', 50, 'Color', color);
    
    % Overlay unit circle once (using a tag to check if it exists)
    if isempty(findobj(gca, 'Tag', 'kp_circle'))
        theta = linspace(0, 2*pi, 200);
        plot(cos(theta), sin(theta), '--k', 'LineWidth', 1.2, 'Tag', 'kp_circle');
        text(1.1, 1.1, '$k_p$', 'Interpreter', 'latex', 'FontSize', 12);
    end
    
    % Local axis formatting
    axis equal square;
    grid on;
    xlabel('$k_x / k_p$', 'Interpreter', 'latex');
    ylabel('$k_y / k_p$', 'Interpreter', 'latex');
end