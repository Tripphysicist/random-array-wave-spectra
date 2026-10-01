function [x, y, l2] = randomParticleLocationsFromL2(numberOfParticles, ensembles, targetL2)
% Generates particle locations from a Gaussian distribution scaled so that
% the mean squared interparticle distance converges to the targetL2.

if nargin < 3
    error('Requires three inputs: numberOfParticles, ensembles, targetL2');
end

% Calculate required standard deviation from the target l2
stdDev = sqrt(targetL2 / 4);

% Generate N by M coordinate matrices
x = stdDev .* randn(numberOfParticles, ensembles);
y = stdDev .* randn(numberOfParticles, ensembles);

% Calculate combinations for all particle pairs
C = nchoosek(1:numberOfParticles, 2);

% Vectorized calculation of squared distances for all pairs across all ensembles
xDist2 = (x(C(:,1), :) - x(C(:,2), :)).^2;
yDist2 = (y(C(:,1), :) - y(C(:,2), :)).^2;

dist = xDist2 + yDist2;

% Mean squared distance across all pairs and ensembles
%l2 = mean(dist, 'all');

% Mean squared distance across all pairs for each ensemble
l2 = mean(dist, 1);


end