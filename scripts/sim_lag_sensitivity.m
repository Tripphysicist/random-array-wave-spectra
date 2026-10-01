% Setting up directional spectra analysis based on WAFO routines
%dMethods = {'MLM','IMLM','EMEM','BDM','EWDM','music','idm','ga','bdm_nnls'};
dMethods = {'BDM'};
numMethods = numel(dMethods);
N = 7; 
close all;
clearvars -except N dMethods numMethods;
tic;

%% Setup global parameters, time vector, and time lag levels
% Lags defined as a percentage of peak period (e.g., 0.5%, 1%, 2.5%, 5%)
lagFractions = [0.001, 0.005, 0.01, 0.05, 0.10]; 
numLags = numel(lagFractions);

waveFieldEnsemble = 30;
particleEnsemble  = 30;
samplingMinutes = 30;
simDuration = samplingMinutes * 60; 
dt = 0.25; 
Nt = (1/dt) * simDuration;
simTime = (0:Nt-1)' * dt; 

%% Prescribe directional spectrum
Hm0 = 1.5;
Tp = 6; % Peak Period
depth = inf;
parametricSpectra  = jonswap(2*pi*linspace(0.04,0.6,64),[Hm0 Tp]);
dTheta  = 90;
meanDirection = 0;
directionalSpread1  = 15;
directionalSpread2  = 15;
dataoptions = [directionalSpread1 directionalSpread2 0.52 5 -2.5 0 1 inf];
parametricDirectionalDistribution = spreading(dTheta, 'cos2', meanDirection,...
    dataoptions, parametricSpectra.w, 1);
parametricDirectionalSpectrum = mkdspec(parametricSpectra,parametricDirectionalDistribution);

%% Generate spatial arrays based on target L2
numberOfParticles = N; 
lambdap = (2*pi)/w2k((2*pi)/Tp, 0, depth);
targetL2 = (lambdap*0.03)^2; 
[x_coords, y_coords, l2] = randomParticleLocationsFromL2(numberOfParticles, particleEnsemble, targetL2);

%% Start Time Lag Loop
for lIdx = 1:numLags
    currLagFrac = lagFractions(lIdx);
    lagMagnitude = currLagFrac * Tp; % Lag in seconds
    fprintf('\n--- Processing Time Lag: %.1f%% of Tp (%.3fs) ---\n', currLagFrac*100, lagMagnitude);

    % Preallocate
    sampledDirectionConfAll = zeros(dTheta, waveFieldEnsemble, particleEnsemble, numMethods, 'single');
    sampledDirectionalSpecConfAll = zeros(dTheta, length(parametricSpectra.w), waveFieldEnsemble, particleEnsemble, numMethods, 'single');

    %% Main ensemble loop
    for ens = 1:waveFieldEnsemble
        disp(['Running ensemble ' num2str(ens)])
        
        [~, currentEnsemblePhases] = simulateSensorTimeseries(parametricDirectionalSpectrum, 0, 0, simTime(1));
        
        for xx = 1:particleEnsemble
            x_pos = x_coords(:,xx)';
            y_pos = y_coords(:,xx)';
            
            % 1. Create a unique time lag for each sensor
            % This simulates asynchronous logging or "clock drift"
            delta_t = (rand(1, N) - 0.5) * 2 * lagMagnitude;
            
            % 2. Generate time series for each sensor with its specific lag
            % We must loop per sensor because simulateSensorTimeseries applies 
            % a single time vector to all positions.
            sensorTimeSeries = zeros(Nt, N);
            for s = 1:N
                % Shift the time vector for THIS sensor
                shiftedTime = simTime + delta_t(s);
                
                tempTS = simulateSensorTimeseries(parametricDirectionalSpectrum, ...
                    x_pos(s) - mean(x_pos), ...
                    y_pos(s) - mean(y_pos), ...
                    shiftedTime, ...
                    currentEnsemblePhases);
                sensorTimeSeries(:,s) = tempTS;
            end
            
            % Use a relative noise floor (1%) instead of an absolute 1mm floor
            sigma_signal = std(sensorTimeSeries(:));
            noise_floor = randn(size(sensorTimeSeries)) * (sigma_signal * 0.01); 

            sampledData = [simTime, sensorTimeSeries + noise_floor];
            
            % 3. Setup WAFO Position Input
            sensorTypes = repmat(sensortypeid('n'), N, 1);
            positionData = [(x_pos - mean(x_pos))', (y_pos - mean(y_pos))', zeros(N,1)];
            positionInput = [positionData, sensorTypes, ones(N, 1)];
            
            % Process through methods
            for mIdx = 1:numMethods
                spec = dat2dspec(sampledData, positionInput, depth, (length(parametricSpectra.w)-1)*2, dTheta, dMethods{mIdx});
                
                sampledDirectionalSpecConfAll(:,:,ens,xx,mIdx) = single(spec.S);
                sampledDirectionConfAll(:,ens,xx,mIdx) = single(trapz(spec.w, spec.S'));
                
                if ens == 1 && xx == 1 && mIdx == 1
                    sampledDirectionalSpectrum_base = spec;
                    sampledDirectionalSpectrum_base.S = []; 
                end
            end
        end
    end

    %% Process and Save
    sDir = ['/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/'...
        'My Drive/Work/analysis/directionalSpectra/data/test/lag/'];
    if ~exist(sDir, 'dir'), mkdir(sDir); end

    for mIdx = 1:numMethods
        currentMethod = dMethods{mIdx};
        sampledDirectionalSpecConf = sampledDirectionalSpecConfAll(:,:,:,:,mIdx);
        sampledDirectionConf = sampledDirectionConfAll(:,:,:,mIdx);
        
        sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
        sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
        
        saveName = sprintf('%s%s_N%02d_lag%.3fTp_ens%02d_cfg%02d.mat', ...
            sDir, currentMethod, N, currLagFrac, waveFieldEnsemble, particleEnsemble);
            
        save(saveName, 'sampledDirectionalSpecConf', 'sampledDirectionConf', ...
                       'sampledEnsembleSpectrum', 'sampledEnsembleDistribution', ...
                       'parametricDirectionalSpectrum', 'N', 'targetL2', ...
                       'x_coords', 'y_coords', 'currLagFrac', 'Tp');
                       
        disp(['Saved: ' saveName]);
    end
end
toc