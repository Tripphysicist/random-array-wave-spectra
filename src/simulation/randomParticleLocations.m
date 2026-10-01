function [x, y, l2] = randomParticleLocations(numberOfParticles,timeSteps,ensembles,scalar);
% function [x, y, l2] = randomParticleLocations(numberOfParticles,timeSteps,ensembles,scalar);
%
% random particle dispersion
%
%

%%


% % find the value of exponential decay of the spectrum, nu, from input
% spectrum

% wavenumberSpectrum = spec2spec(inputSpectrum, 'k1d');

% [~, findex] = max(wavenumberSpectrum.S); %find the peak of the spectrum
% lambdap = (2*pi)./wavenumberSpectrum.k(findex); %peak wavelength is the length scale

% % fit to the tail (2 x kp) in log space
% pfit = polyfit(log(wavenumberSpectra.w(findex*2:end)),log(wavenumberSpectra.S(findex*2:end)),1);
% nu = abs(pfit(1)); %the absolute value of the slope is nu

% the second coefficient comes from the dispersion relation which can be
% measured or assumed
% nu = 3; % for Phillips spectrum f^-5
% %nu = 5/2; % for Zakarov spectrum f^-4
% alpha = 1/2;
% lambda = 4 ./(2+2-nu);
% mu = lambda - 1./alpha;

% each time step, update the x and y coordinates

if nargin < 3
    error('need more inputs')
end

C = nchoosek(1:numberOfParticles, 2); % all possible combinations of particles

%clr = lansey(numberOfParticles);

%
% want to generate and evolve the movement of particles. the pdf of the
% motion can be Gaussian (random) or based on another pdf of choice.
%
%
% based on Cartesian grid x, y
% initialize location

x = zeros(numberOfParticles,timeSteps,ensembles);
y = zeros(numberOfParticles,timeSteps,ensembles);
%d = zeros(timeSteps,ensembles);
dist = zeros(length(C),timeSteps,ensembles);

% figure
% subplot(2,1,1)
% plot(-100:100,zeros(size(-100:100)),'k--')
% hold on
% plot(zeros(size(-100:100)),-100:100,'k--')
% grid on


for ee = 1:ensembles
    for i = 1:timeSteps-1
        for jj = 1:numberOfParticles

            x(jj,i+1,ee) = x(jj,i,ee) + randn.*scalar;
            y(jj,i+1,ee) = y(jj,i,ee) + randn.*scalar;
            
            %track the distance between particle 1 and 2 ad 1 and 3
            %there are equation for conbinations of 2 particles is 
            % C = nchoosek(numberOfParticles, 2)
            % the mean square inter-particle distance (diameter of particle cluster)
            for nn = 1:length(C)
                xDist2 = (x(C(nn,1)+1,i+1,ee) - x(C(nn,2),i+1,ee)).^2;
                yDist2 = (y(C(nn,1)+1,i+1,ee) - y(C(nn,2),i+1,ee)).^2;
                dist(nn,i+1,ee) = xDist2 + yDist2; 
            end
            %         plot(x(jj),y(jj),'.','Color',clr(jj,:),'markerSize',25)
            %         shg
            %         hold off
            %         plot(-100:100,zeros(size(-100:100)),'k--')
            %         hold on
            %         plot(zeros(size(-100:100)),-100:100,'k--')
            %         grid on
        end
    end
end
l2 = mean(mean(dist,3));

% %% plot mean square distance
% plot(l2,'LineWidth',3)
% legend('$\langle l^2(t) \rangle$','Interpreter','latex')
% fontsize(20,20,20,20)
% xlabel('time')
% ylabel('distance')
% grid on
% 
%
end