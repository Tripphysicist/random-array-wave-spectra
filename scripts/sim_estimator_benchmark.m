% Setting up directional spectra analysis based on WAFO routines
dMethods = {'MLM','IMLM','EMEM','BDM','EWDM','music','idm','ga','bdm_nnls'};

%dMethods = {'MLM','IMLM'};
numMethods = numel(dMethods);
N = 15; % Explicitly define sensor count here

close all;
clearvars -except N dMethods numMethods;
tic;

%% Setup global parameters and time vector
waveFieldEnsemble = 30;
particleEnsemble  = 30;

samplingMinutes = 30;
simDuration = samplingMinutes * 60; % [seconds]
dt = 0.25; % sampling frequency (seconds)
Nt = (1/dt) * simDuration;
simTime = (0:Nt-1)' * dt; % Column vector of time steps

%% Prescribe directional spectrum
Hm0 = 1.5;
Tp = 6;
lambdap = (2*pi)/w2k((2*pi)/Tp);
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

%% Estimator setup - Add 4th dimension for 'method'
% [theta x ensemble x config x method]
sampledDirectionConfAll = zeros(dTheta, waveFieldEnsemble, particleEnsemble, numMethods, 'single');
% [theta x freq x ensemble x config x method]
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
        
        % Simulate exact time series ONCE per configuration
        sensorTimeSeries = simulateSensorTimeseries(parametricDirectionalSpectrum, ...
            sensorLocationXdemean, ...
            sensorLocationYdemean, ...
            simTime, ...
            currentEnsemblePhases);
            
        sampledData = [simTime, sensorTimeSeries];
        sensorTypes = repmat(sensortypeid('n'), numSensors, 1);
        positionData   = [sensorLocationXdemean', sensorLocationYdemean', zeros(numSensors,1)];
        depth = inf;
        bfuncs   = ones(numSensors, 1);
        positionInput = [positionData sensorTypes, bfuncs];
        
        % Process the identical time series through all 6 methods
        for mIdx = 1:numMethods
            spec = dat2dspec(sampledData, positionInput, depth, (length(parametricSpectra.w)-1)*2, dTheta, dMethods{mIdx});
            
            sampledDirectionalSpecConfAll(:,:,ens,xx,mIdx) = single(spec.S);
            sampledDirectionConfAll(:,ens,xx,mIdx) = single(trapz(spec.w, spec.S'));
            
            % Capture the theta/w vectors from the struct on the first pass
            if ens == 1 && xx == 1 && mIdx == 1
                sampledDirectionalSpectrum_base = spec;
                sampledDirectionalSpectrum_base.S = []; % Clear massive matrix to save overhead
            end
        end
    end
end

%% Process and Save Files per Method
parametricDistribution = trapz(parametricDirectionalSpectrum.w, parametricDirectionalSpectrum.S');
sDir = '/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/My Drive/Work/analysis/directionalSpectra/data/test/methods/N15/';

for mIdx = 1:numMethods
    currentMethod = dMethods{mIdx};
    
    % Extract the specific 3D data for the current method
    sampledDirectionalSpecConf = sampledDirectionalSpecConfAll(:,:,:,:,mIdx);
    sampledDirectionConf = sampledDirectionConfAll(:,:,:,mIdx);
    
    sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
    sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
    
    % Reconstruct the WAFO output struct for downstream plotting compatibility
    sampledDirectionalSpectrum = sampledDirectionalSpectrum_base;
    
    saveName = sprintf('%s%s_sensors%02d__targetL2_%0.4f_ensembles%02d_configs%02d.mat', ...
        sDir, currentMethod, N, targetL2, waveFieldEnsemble, particleEnsemble);
        
    save(saveName, 'sampledDirectionalSpecConf', 'sampledDirectionConf', ...
                   'sampledEnsembleSpectrum', 'sampledEnsembleDistribution', ...
                   'parametricDistribution', 'parametricDirectionalSpectrum', ...
                   'sampledDirectionalSpectrum', 'N', 'targetL2', 'x', 'y', 'l2');
                   
    disp(['Saved: ' saveName]);
end
toc