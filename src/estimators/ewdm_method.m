function [S2D, w, theta] = ewdm_method(xn, pos, nfft, nt)
    % 1. Prepare data for export
    t = xn(:, 1);
    Fs = 1/mean(diff(t));
    signal = xn(:, 2:end);
    x = pos(:, 1);
    y = pos(:, 2);
    
    save('tmp_ewdm_input.mat', 't', 'signal', 'x', 'y', 'nfft', 'nt');
    
    % 2. Call Python Bridge
    pythonPath = '/Users/tripp/anaconda3/envs/py310_intel/bin/python';
    scriptPath = '/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/My Drive/Work/analysis/directionalSpectra/methods and 8-m array/extended-wdm/ewdm/ewdm_bridge.py';
    
    cmd = sprintf('%s "%s"', pythonPath, scriptPath);
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('EWDM Bridge Failed: %s', cmdout);
    end
    
    % 3. Load results
    res = load('tmp_ewdm_output.mat');
    S_log = res.S; % [frequency_hz x direction_deg]
    f_hz_log = res.f_hz(:);
    theta_deg = res.theta_deg(:);
    
    % --- Interpolation to Linear Grid ---
    f_linear = (0:nfft/2)' * Fs / nfft;
    S_linear = zeros(length(f_linear), length(theta_deg));
    for i = 1:length(theta_deg)
        S_linear(:,i) = interp1(f_hz_log, S_log(:,i), f_linear, 'pchip', 0);
    end
    
    % Convert Hz vector to rad/s for WAFO's structural requirement
    w = f_linear * (2 * pi);
    
    % 1. Direction: Meteorological (Coming From) -> Mathematical (Going To)
    theta_deg_math = theta_deg + 180;
    
    % Convert to radians and wrap perfectly to [-pi, pi]
    theta = theta_deg_math * (pi / 180);
    theta = atan2(sin(theta), cos(theta));
    
    % Sort angles monotonically for WAFO's integration routines
    [theta, sortIdx] = sort(theta);
    
    % Format S2D (Transpose to [theta x w])
    % NO SCALING. The data is perfectly aligned.
    S2D = S_linear(:, sortIdx)';
    
    % Cleanup
    delete('tmp_ewdm_input.mat');
    delete('tmp_ewdm_output.mat');
end