% Setting up directional spectra analysis based on WAFO routines


scale = logspace(-4,1,11);
for bigL = [1:3 10:11]
    close all
    clearvars -except scale bigL
    tic

    %% Setup global parameters and time vector

    waveFieldEnsemble = 30;
    particleEnsemble = 100;
    %targetL2 = 0.2903; % Target root-mean-square distance parameter

    N=10; %number of sensors

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

    %% Define this once before your ensemble loop begins
    samplingMinutes = 30;
    simDuration = samplingMinutes * 60; % [seconds]
    dt = 0.25; % sampling frequency (seconds)
    Nt = (1/dt) * simDuration;
    simTime = (0:Nt-1)' * dt; % Column vector of time steps


    %% Estimator setup
    dMethods = {'BDM'};
    sampledDirectionConf = zeros(dTheta, waveFieldEnsemble, particleEnsemble);
    sampledDirectionalSpecConf = zeros(dTheta, length(parametricSpectra.w), waveFieldEnsemble, particleEnsemble);

    %% Generate ALL random spatial arrays based on target L2

    targetL2 = (lambdap*scale(bigL))^2; % Target root-mean-square distance parameter

    % Using the updated random generator function directly
    [x, y, l2] = randomParticleLocationsFromL2(N, particleEnsemble, targetL2);

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

    sDir = ['/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/' ...
        'My Drive/Work/analysis/directionalSpectra/data/test/targetL2/N10/'];

    % Use %02d for zero-padding the sensor count to ensure correct directory sorting
    saveName = sprintf('%ssensors%02d_scale_%02d_ensembles%02d_configs%02d.mat', ...
        sDir, N, bigL, waveFieldEnsemble, particleEnsemble);
    save(saveName)
    toc
end


%% plot from the different data files
%
clear
close all
dataFolder = '/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/My Drive/Work/analysis/directionalSpectra/data/test/targetL2adapt/';
datadir = dir(dataFolder);
datadir=datadir(3:end);

%
for loopi=1:numel(datadir)
    load([dataFolder  datadir(loopi).name])
    %ensemble average
    sampledEnsembleSpectrum = squeeze(mean(sampledDirectionalSpecConf,3));
    sampledEnsembleDistribution = squeeze(mean(sampledDirectionConf,2));
    parametricDistribution=trapz(parametricDirectionalSpectrum.w,parametricDirectionalSpectrum.S');
    %plotting
    %ll2(loopi) = l2(end);
    figure(1)
    subplot(3,4,loopi)
    plot(sampledDirectionalSpectrum.theta,sampledEnsembleDistribution)
    hold on
    plot(parametricDirectionalSpectrum.theta,parametricDistribution,'k--','LineWidth',4)
    plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'b','LineWidth',2)
    grid on
    xlabel('direction [deg]')
    ylabel('S(\theta) [m^z/deg]')
%    fontsize(18,18,18,18)
    hold on

    [~, maxIndex]=max(sampledEnsembleDistribution);
    %pg22=sum(abs(parametricDirectionalSpectrum.theta(maxIndex))>pi/16); %percent outside +-11.25 deg
    stddir(loopi) = std(sampledDirectionalSpectrum.theta(maxIndex))*(180/pi);
    %    f = kde(parametricDirectionalSpectrum.theta(maxIndex));
    %    figure(2)
    %    plot(cell2mat(f.x),f.f,'LineWidth',3); hold on
    %    histogram(parametricDirectionalSpectrum.theta(maxIndex)*(180/pi),11); hold on

    %mean energy direction
    sumsumerror(loopi) = sum(sum(abs(sampledEnsembleDistribution-parametricDistribution')));
%    rmsd(loopi)=l2(l2index);


    meanDS1 = mean(sampledDirectionalSpecConf,3);
    meanDS2 = mean(meanDS1,4);


    figure(2)
    subplot(3,4,loopi)
    pcolor(sampledDirectionalSpectrum.theta,sampledDirectionalSpectrum.w,meanDS2')
    shading("interp")
%    colormap(lansey);

    xlabel('direction [deg]')
    ylabel('\omega [Hz]')

%    fontsize(18,18,18,18)
    hold on

end
%
figure(3)
subplot(2,1,1)
semilogx(scale,stddir,'o-','LineWidth',3)
grid on
xlabel('l_{rms}/\lambda_p','Interpreter','tex')
ylabel('\sigma(\theta_p) [\circ]')
%fontsize(18,18,18,18)

subplot(2,1,2)
semilogx(scale,sumsumerror,'o-','LineWidth',3)
grid on
xlabel('l_{rms}/\lambda_p','Interpreter','tex')
ylabel('sum error [m^2]')
%fontsize(18,18,18,18)

%%

%% Two metrics, bias in peak direction and mean (maybe sum?) energy
clear
close all
dataFolder = ['/Users/tripp/Library/CloudStorage/GoogleDrive-tripphysicist@gmail.com/...'...
    'My Drive/Work/analysis/directionalSpectra/data/test/targetL2/'];
datadir = dir(dataFolder);
datadir=datadir(3:end);

count = 0;
for loopi=1:numel(datadir)
    load([dataFolder  datadir(loopi).name])

    count = count + 1;

    %location of peak direction
    [~, maxIndex]=max(sampledEnsembleDistribution);
    pg22=sum(abs(sampledDirectionalSpectrum.theta(maxIndex))>pi/16); %percent outside +-11.25 deg
    stddir = std(sampledDirectionalSpectrum.theta(maxIndex))*(180/pi);

    %mean energy direction
    sumerror = sum(abs(sampledEnsembleDistribution-parametricDistribution'));

    % plot polygons from arrays
    %[B,I]=sort(abs(sampledDirectionalSpectrum.theta(maxIndex)));
    [B,I]= sort(sumerror);

    meanB(loopi)=mean(B);
    stdB(loopi)=std(B);

    clr=colororder;
    clr2 = lansey(10);
    sc = 3;
    %for i=1:3 % top 3
    for i=[1 30] % bottom 3
        % Find centroid
%         xp = x(:,l2index,I(i)); older file
%         yp = y(:,l2index,I(i));
        xp = x(:,I(i));
        yp = y(:,I(i));
        centroidX = mean(xp);
        centroidY = mean(yp);

        % Calculate angles
        angles = atan2(yp - centroidY, xp - centroidX);

        % Sort points by angles
        [sortedAngles, sortIndices] = sort(angles);

        % Reorder points
        reorderedX = xp(sortIndices);
        reorderedY = yp(sortIndices);

        % Create and plot polygon
        polygon = polyshape([reorderedX, reorderedY]);
        [xv, yv] = boundary(polygon);


        figure(1)

        subplot(2,3,count)
        pp = plot(polygon);

        if i==1
            pp.FaceColor = clr(1,:);
            pp.EdgeColor = clr(1,:);
            hold on;
            plot(centroidX, centroidY, 'r+', 'MarkerSize', 10,'Color',clr(1,:)); % Plot the centroid
            plot(xv, yv, 'Marker' ,'.', 'MarkerSize', 20,'Color',clr(1,:))
        else
            pp.FaceColor = clr(2,:);
            pp.EdgeColor = clr(2,:);
            hold on;
            plot(centroidX, centroidY, 'r+', 'MarkerSize', 10,'Color',clr(2,:)); % Plot the centroid
            plot(xv, yv, 'Marker' ,'.', 'MarkerSize', 20,'Color',clr(2,:))
        end
        axis([-sc sc -sc sc]); % Ensure equal aspect ratio for a visually correct polygon
        grid on
        title(num2str(i))
        xlabel('x [m]')
        ylabel('y [m]')
        %hold off;



        %     subplot(1,3,2)
        %     coarrayAnalysisPlots(xp,yp,Tp,i)
        %
        %     subplot(1,3,3)
        %     if count==1
        %         plot(parametricDirectionalSpectrum.theta,parametricDistribution,'k--','LineWidth',4)
        %     end
        %     hold on
        %     plot(sampledDirectionalSpectrum.theta,sampledEnsembleDistribution(:,I(i)),'LineWidth',3,'Color',clr(count,:))
        %     %hold off;
        grid on
        %     legend
        %     xlabel('\theta [\circ]')
        %     ylabel('E(\theta) [m^2/deg]')
        fontsize(18,18,18,18)
        shg

    end
    figure(3)
    subplot(2,3,count)
    shadedErrorBar(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),std(sampledEnsembleDistribution'),'lineProps','b')
    %plot(sampledDirectionalSpectrum.theta,sampledEnsembleDistribution)
    hold on
    plot(parametricDirectionalSpectrum.theta,parametricDistribution,'k--','LineWidth',4)
    %plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'b','LineWidth',2)
    text(-3,0.06,['N = ' num2str(N(loopi))])
    grid on
    xlabel('direction [deg]')
    ylabel('S(\theta) [m^z/deg]')
    fontsize(18,18,18,18)
    axis([-pi pi 0 0.08])
%     subplot(2,3,11)
%     if loopi==1
%         plot(parametricDirectionalSpectrum.theta,parametricDistribution,'k--','LineWidth',4)
%         hold on
%         plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'Color',clr2(loopi,:),'LineWidth',2)
%     else
%         plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'Color',clr2(loopi,:),'LineWidth',2)
%     end
%     axis([-pi pi 0 0.08])
%     grid on
%     xlabel('direction [deg]')
%     ylabel('S(\theta) [m^z/deg]')
%     fontsize(18,18,18,18)

    figure(2)
    if loopi==1
        plot(parametricDirectionalSpectrum.theta,parametricDistribution,'k--','LineWidth',4)
        hold on
        plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'Color',clr2(loopi,:),'LineWidth',2)
    else
        plot(sampledDirectionalSpectrum.theta,mean(sampledEnsembleDistribution'),'Color',clr2(loopi,:),'LineWidth',2)
    end
    grid on
    xlabel('direction [deg]')
    ylabel('S(\theta) [m^z/deg]')
    fontsize(18,18,18,18)
end

figure
shadedErrorBar(scale,meanB,stdB)