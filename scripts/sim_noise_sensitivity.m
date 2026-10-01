% Setting up directional spectra analysis based on WAFO routines
%dMethods = {'MLM','IMLM','EMEM','BDM','EWDM'};%,'music','idm','ga','bdm_nnls'};
dMethods = {'BDM'};
numMethods = numel(dMethods);
N = 7; % Explicitly define sensor count here
close all;
clearvars -except N dMethods numMethods;
tic;

%% Setup global parameters, time vector, and noise floors
% noiseLevels_meters = [0.005, 0.001, 0.05, 0.10, 0.5, 1.0]; % Physical floors (m)
noiseLevels_meters = 0.50;
numNoise = numel(noiseLevels_meters);
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
lambdap = (2*pi)/w2k((2*pi)/Tp, 0, depth);
parametricSpectra  = jonswap(2*pi*linspace(0.04,0.6,64),[Hm0 Tp]);
dTheta  = 90;
meanDirection = 0;
directionalSpread1  = 15;
directionalSpread2  = 15;
dataoptions = [directionalSpread1 directionalSpread2 0.52 5 -2.5 0 1 inf];
parametricDirectionalDistribution = spreading(dTheta, 'cos2', meanDirection,...
    dataoptions, parametricSpectra.w, 1);
parametricDirectionalSpectrum = mkdspec(parametricSpectra,parametricDirectionalDistribution);

%% Generate ALL random spatial arrays based on target L2
numberOfParticles = N; 
targetL2 = (lambdap*0.03)^2; 
[x, y, l2] = randomParticleLocationsFromL2(numberOfParticles, particleEnsemble, targetL2);

%% Start Noise Loop
for nIdx = 1:numNoise
    currSigma = noiseLevels_meters(nIdx);
    fprintf('\n--- Processing Noise Level: %.3f m ---\n', currSigma);

    % Preallocate for this specific noise level
    sampledDirectionConfAll = zeros(dTheta, waveFieldEnsemble, particleEnsemble, numMethods, 'single');
    sampledDirectionalSpecConfAll = zeros(dTheta, length(parametricSpectra.w), waveFieldEnsemble, particleEnsemble, numMethods, 'single');

    %% Main ensemble loop
    for ens = 1:waveFieldEnsemble
        disp(['Running ensemble ' num2str(ens)])
        
        % Generate random wave phases ONCE per ensemble
        [~, currentEnsemblePhases] = simulateSensorTimeseries(parametricDirectionalSpectrum, 0, 0, simTime(1));
        
        for xx = 1:particleEnsemble
            sensorLocationX = x(:,xx)';
            sensorLocationY = y(:,xx)';
            numSensors = numel(sensorLocationX);
            
            sensorLocationXdemean = sensorLocationX - mean(sensorLocationX);
            sensorLocationYdemean = sensorLocationY - mean(sensorLocationY);
            
            % 1. Simulate exact time series
            sensorTimeSeries = simulateSensorTimeseries(parametricDirectionalSpectrum, ...
                sensorLocationXdemean, ...
                sensorLocationYdemean, ...
                simTime, ...
                currentEnsemblePhases);
            
            % 2. Add random noise for this realization
            noisyTimeSeries = sensorTimeSeries + randn(size(sensorTimeSeries)) * currSigma;
            sampledData = [simTime, noisyTimeSeries];
            
            % Setup Position Input
            sensorTypes = repmat(sensortypeid('n'), numSensors, 1);
            positionData   = [sensorLocationXdemean', sensorLocationYdemean', zeros(numSensors,1)];
            bfuncs   = ones(numSensors, 1);
            positionInput = [positionData sensorTypes, bfuncs];
            
            % Process through all methods
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

    %% Process and Save Files for this Noise Level
    parametricDistribution = trapz(parametricDirectionalSpectrum.w, parametricDirectionalSpectrum.S');
    sDir = '/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/My Drive/Work/analysis/directionalSpectra/data/test/noise/';
    if ~exist(sDir, 'dir'), mkdir(sDir); end

    for mIdx = 1:numMethods
        currentMethod = dMethods{mIdx};
        
        sampledDirectionalSpecConf = sampledDirectionalSpecConfAll(:,:,:,:,mIdx);
        sampledDirectionConf = sampledDirectionConfAll(:,:,:,mIdx);
        
        sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
        sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
        
        sampledDirectionalSpectrum = sampledDirectionalSpectrum_base;
        
        % Filename now includes noise level for clarity
        saveName = sprintf('%s%s_N%02d_noise%.3f_L2_%0.4f_ens%02d_cfg%02d.mat', ...
            sDir, currentMethod, N, currSigma, targetL2, waveFieldEnsemble, particleEnsemble);
            
        save(saveName, 'sampledDirectionalSpecConf', 'sampledDirectionConf', ...
                       'sampledEnsembleSpectrum', 'sampledEnsembleDistribution', ...
                       'parametricDistribution', 'parametricDirectionalSpectrum', ...
                       'sampledDirectionalSpectrum', 'N', 'targetL2', 'x', 'y', 'l2', 'currSigma');
                       
        disp(['Saved: ' saveName]);
    end
end
toc