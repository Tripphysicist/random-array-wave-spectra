% Setting up directional spectra analysis based on WAFO routines

%N=[3,4,5,6,7,8,9,10,15,20];
N=[10];

for bigL = 1:numel(N)
    close all
    clearvars -except N bigL
    tic

    %% Setup global parameters and time vector
    
    waveFieldEnsemble = 30;
    particleEnsemble = 30;
    %targetL2 = 0.2903; % Target root-mean-square distance parameter

    %% Define this once before your ensemble loop begins
    samplingMinutes = 30;
    simDuration = samplingMinutes * 60; % [seconds]
    dt = 0.25; % sampling frequency (seconds)
    %dt = 0.10;
    Nt = (1/dt) * simDuration;
    simTime = (0:Nt-1)' * dt; % Column vector of time steps

        %% Prescribe directional spectrum
        Hm0 = 1.5;
        Tp = 6;
        wp = (2*pi)/Tp;
        lambdap = (2*pi)/w2k(wp);
        parametricSpectra  = jonswap(2*pi*linspace(0.04,0.6,64),[Hm0 Tp]);
    
        dTheta  = 90;
        meanDirection = 0;
        directionalSpread1  = 15;
        directionalSpread2  = 15;
        dataoptions = [directionalSpread1 directionalSpread2 wp 5 -2.5 0 1 inf];
    
    %     parametricDirectionalDistribution = spreading(dTheta, 'bimodal', meanDirection,...
    %         dataoptions, parametricSpectra.w, 1);
    
        parametricDirectionalDistribution = spreading(dTheta, 'cos2', meanDirection,...
            dataoptions, parametricSpectra.w, 1);
    
        parametricDirectionalSpectrum = mkdspec(parametricSpectra,parametricDirectionalDistribution);

    %% Estimator setup
    dMethods = {'BDM'};
    sampledDirectionConf = zeros(dTheta, waveFieldEnsemble, particleEnsemble);
    sampledDirectionalSpecConf = zeros(dTheta, length(parametricSpectra.w), waveFieldEnsemble, particleEnsemble);

    %% Generate ALL random spatial arrays based on target L2
    numberOfParticles = N(bigL);
    
    targetL2 = (lambdap*0.01)^2; % Target root-mean-square distance parameter

    % Using the updated random generator function directly
    [x, y, l2] = randomParticleLocationsFromL2(numberOfParticles, particleEnsemble, targetL2);

    %% Main ensemble loop
    for ens = 1:waveFieldEnsemble
        disp(['running ensemble ' num2str(ens)])

        % Generate random wave phases ONCE per ensemble
        % We call the simulator with empty coords just to extract the Phi_flat array
        [~, currentEnsemblePhases] = simulateSensorTimeseries(parametricDirectionalSpectrum, 0, 0, simTime(1));

        for xx = 1:particleEnsemble
            % Extract sensor locations for this configuration
            sensorLocationX = x(:,xx)';
            sensorLocationY = y(:,xx)';
            numSensors = numel(sensorLocationX);

            % Demean to preserve floating point precision
            sensorLocationXdemean = sensorLocationX - mean(sensorLocationX);
            sensorLocationYdemean = sensorLocationY - mean(sensorLocationY);
            % Simulate exact time series using the fixed ensemble phases
            sensorTimeSeries = simulateSensorTimeseries(parametricDirectionalSpectrum, ...
                sensorLocationXdemean, ...
                sensorLocationYdemean, ...
                simTime, ...
                currentEnsemblePhases);
            %% calculate spectrum from time series
%             for i = 1:numSensors
%                 data = [simTime, sensorTimeSeries(:,i)];
%                 sampledSpectrumTemp = dat2spec(data,200);
%                 if i == 1
%                     sampledSpectrum = sampledSpectrumTemp;
%                 end
%                 sensorNo = num2str(i);
%                 sampledSpectrum.(['S' sensorNo]) = sampledSpectrumTemp.S;
%                 sampledSpectrum.bigS(:,i) = sampledSpectrumTemp.S;
%             end
%             sampledSpectrum = rmfield(sampledSpectrum, 'S');

            %% calculate directional spectrum from time series
            sampledData = [simTime, sensorTimeSeries];
            sensorTypes = repmat(sensortypeid('n'), numSensors, 1);
            positionData   = [sensorLocationXdemean', sensorLocationYdemean', zeros(numSensors,1)];
            depth = inf;
            bfuncs   = ones(numSensors, 1);
            positionInput = [positionData sensorTypes, bfuncs];
            sampledDirectionalSpectrum = dat2dspec(sampledData, positionInput, depth, (length(parametricSpectra.w)-1)*2, dTheta, dMethods{1});
            sampledDirectionalSpecConf(:,:,ens,xx) = sampledDirectionalSpectrum.S;
            sampledDirectionConf(:,ens,xx) = trapz(sampledDirectionalSpectrum.w, sampledDirectionalSpectrum.S');
        end
    end

    %% average over wave field ensembles and save
    sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
    sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
    parametricDistribution=trapz(parametricDirectionalSpectrum.w,parametricDirectionalSpectrum.S');

    % sDir = ['/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/'...
    %    'My Drive/Work/analysis/directionalSpectra/data/test/targetL2/N10/'];

% Use %02d for zero-padding the sensor count to ensure correct directory sorting
    saveName = sprintf('%ssensors%02d__targetL2_%0.4f_ensembles%02d_configs%02d.mat', ...
    sDir, N(bigL), targetL2, waveFieldEnsemble, particleEnsemble);
    save(saveName)
    toc
end

%%

%% Extract and Plot Sampled Spreading Function
% 1. Average over ensembles (dim 3) and configurations (dim 4)
% Resulting matrix is [90 x 64] (theta x w)         
S_mean_2D = squeeze(mean(mean(sampledDirectionalSpecConf, 3), 4));

% 2. Calculate the 1D frequency spectrum S(w) by integrating over theta (dim 1)
% Ensure your 'theta' vector is in radians for trapz
S_w = trapz(parametricDirectionalSpectrum.theta, S_mean_2D, 1); % Result is [1 x 64]

% 3. Divide out the energy to get D(theta, w)
D_sampled = S_mean_2D ./ S_w;

% 4. Noise Safeguard: Zero out D(theta, w) where there is no physical energy
% (e.g., less than 0.5% of the peak energy) to prevent division-by-near-zero static
% energy_threshold = 0.00001 * max(S_w);
% valid_freqs = S_w > energy_threshold;
% D_sampled(:, ~valid_freqs) = 0;

% 5. Plot the result

%figure('Name', 'Bimodal Spreading Surface', 'Color', 'w', 'Position', [100 100 800 600]);
figure
% Extract vectors and matrix
w_rad = parametricDirectionalSpectrum.w;
theta_deg = parametricDirectionalSpectrum.theta * (180/pi);
D_matrix = parametricDirectionalDistribution.S';

% Option A: Smooth surface plot
subplot(1,2,1)
surf(theta_deg, w_rad, D_matrix, 'EdgeColor', 'none');
view(0, 90); % Top-down view
colormap(turbo); % Turbo provides much better contrast for lobes than Jet
ylim([0.4 3.5])

% Format axes
xlabel('Direction \theta [deg]', 'FontSize', 14);
ylabel('Frequency \omega [rad/s]', 'FontSize', 14);
title('Bimodal Spreading Function D(\theta, \omega)', 'FontSize', 16);
axis tight;
colorbar;

% Add a horizontal line at the peak frequency to mark the split point
yline(wp, 'w--', 'LineWidth', 2, 'Label', '\omega_p');

%figure('Name', 'Sampled Bimodal Spreading', 'Color', 'w', 'Position', [100 100 800 600]);

% Convert axes for plotting (assuming theta and w vectors exist in workspace)
theta_deg = sampledDirectionalSpectrum.theta * (180/pi);
w_rad = sampledDirectionalSpectrum.w;

% Transpose D_sampled to [w x theta] to match the surf orientation
subplot(1,2,2)
surf(theta_deg, w_rad, D_sampled', 'EdgeColor', 'none');
view(0, 90); % Top-down view
colormap(turbo);

% Format axes to match the theoretical plot
xlabel('Direction \theta [deg]', 'FontSize', 14);
ylabel('Frequency \omega [rad/s]', 'FontSize', 14);
title(['Sampled Spreading Function (Method: ', 'BDM', ')'], 'FontSize', 16);
axis tight;
colorbar;
ylim([0.4 3.5])
% Add the peak frequency reference line
wp = (2*pi) / Tp;
yline(wp, 'w--', 'LineWidth', 2, 'Label', '\omega_p');