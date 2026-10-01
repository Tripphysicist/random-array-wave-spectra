% Setting up directional spectra analysis based on WAFO routines
%dMethods = {'MLM','IMLM','EMEM','BDM','music','idm','ga','bdm_nnls'};
dMethods = {'BDM'};
numMethods = numel(dMethods);
N = 7; 
close all;
clearvars -except N dMethods numMethods;
tic;

%% Setup global parameters, time vector, and spatial offset levels
% Offsets defined as a percentage of peak wavelength (e.g., 0.5%, 1%, 2.5%)
offsetFractions = [0.001, 0.005, 0.01, 0.05, 0.1]; 
numOffsets = numel(offsetFractions);

waveFieldEnsemble = 30;
particleEnsemble  = 30;
samplingMinutes = 30;
simDuration = samplingMinutes * 60; 
dt = 0.25; 
Nt = (1/dt) * simDuration;
simTime = (0:Nt-1)' * dt; 

%% Prescribe directional spectrum
Hm0 = 1.5;
Tp = 6;
depth = inf;
kp = w2k((2*pi)/Tp, 0, depth, gravity);
lambdap = (2*pi)/kp; % Peak wavelength

parametricSpectra  = jonswap(2*pi*linspace(0.04,0.6,64),[Hm0 Tp]);
dTheta  = 90;
meanDirection = 0;
directionalSpread1  = 15;
directionalSpread2  = 15;
dataoptions = [directionalSpread1 directionalSpread2 0.52 5 -2.5 0 1 inf];
parametricDirectionalDistribution = spreading(dTheta, 'cos2', meanDirection,...
    dataoptions, parametricSpectra.w, 1);
parametricDirectionalSpectrum = mkdspec(parametricSpectra,parametricDirectionalDistribution);

%% Generate Reported spatial arrays based on target L2
numberOfParticles = N; 
targetL2 = (lambdap*0.03)^2; 
[x_reported, y_reported, l2] = randomParticleLocationsFromL2(numberOfParticles, particleEnsemble, targetL2);

%% Start Spatial Offset Loop
for oIdx = 1:numOffsets
    currOffsetFrac = offsetFractions(oIdx);
    offsetMagnitude = currOffsetFrac * lambdap; % Offset in meters
    fprintf('\n--- Processing Offset Level: %.1f%% of Lp (%.3fm) ---\n', currOffsetFrac*100, offsetMagnitude);

    % Preallocate
    sampledDirectionConfAll = zeros(dTheta, waveFieldEnsemble, particleEnsemble, numMethods, 'single');
    sampledDirectionalSpecConfAll = zeros(dTheta, length(parametricSpectra.w), waveFieldEnsemble, particleEnsemble, numMethods, 'single');

    %% Main ensemble loop
    for ens = 1:waveFieldEnsemble
        disp(['Running ensemble ' num2str(ens)])
        
        [~, currentEnsemblePhases] = simulateSensorTimeseries(parametricDirectionalSpectrum, 0, 0, simTime(1));
        
        for xx = 1:particleEnsemble
            % 1. Get Reported Positions
            x_rep = x_reported(:,xx)';
            y_rep = y_reported(:,xx)';
            
            % 2. Create the "True" (Corrupted) Positions
            % Each sensor is nudged by a random distance up to the offsetMagnitude
            dx = (rand(1, N) - 0.5) * 2 * offsetMagnitude;
            dy = (rand(1, N) - 0.5) * 2 * offsetMagnitude;
            
            x_true = x_rep + dx;
            y_true = y_rep + dy;
            
            % 3. Simulate waves at the TRUE (offset) locations
            sensorTimeSeries = simulateSensorTimeseries(parametricDirectionalSpectrum, ...
                x_true - mean(x_true), ...
                y_true - mean(y_true), ...
                simTime, ...
                currentEnsemblePhases);
            
            % Use a relative noise floor (1%) instead of an absolute 1mm floor
            sigma_signal = std(sensorTimeSeries(:));
            noise_floor = randn(size(sensorTimeSeries)) * (sigma_signal * 0.01); 

            sampledData = [simTime, sensorTimeSeries + noise_floor];
            
            % 4. Build WAFO Position Input using REPORTED positions
            % This simulates the error where the scientist thinks sensors are at (x_rep, y_rep)
            x_rep_demean = x_rep - mean(x_rep);
            y_rep_demean = y_rep - mean(y_rep);
%             x_true_demean = x_true - mean(x_true);
%             y_true_deaman = y_true - mean(y_true);

            sensorTypes = repmat(sensortypeid('n'), N, 1);
            positionData = [x_rep_demean', y_rep_demean', zeros(N,1)];
%             positionDataTrue = [x_true_demean', y_true_deaman', zeros(N,1)];
            positionInput = [positionData, sensorTypes, ones(N, 1)];
%             positionInputTrue = [positionDataTrue, sensorTypes, ones(N, 1)];
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
    sDir = '/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/My Drive/Work/analysis/directionalSpectra/data/test/offset/';
    if ~exist(sDir, 'dir'), mkdir(sDir); end

    for mIdx = 1:numMethods
        currentMethod = dMethods{mIdx};
        sampledDirectionalSpecConf = sampledDirectionalSpecConfAll(:,:,:,:,mIdx);
        sampledDirectionConf = sampledDirectionConfAll(:,:,:,mIdx);
        
        sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
        sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
        
        saveName = sprintf('%s%s_N%02d_offset%.3fLp_L2_%0.4f_ens%02d_cfg%02d.mat', ...
            sDir, currentMethod, N, currOffsetFrac, targetL2, waveFieldEnsemble, particleEnsemble);
            
        save(saveName, 'sampledDirectionalSpecConf', 'sampledDirectionConf', ...
                       'sampledEnsembleSpectrum', 'sampledEnsembleDistribution', ...
                       'parametricDirectionalSpectrum', 'N', 'targetL2', ...
                       'x_reported', 'y_reported', 'currOffsetFrac', 'lambdap');
                       
        disp(['Saved: ' saveName]);
    end
end
toc