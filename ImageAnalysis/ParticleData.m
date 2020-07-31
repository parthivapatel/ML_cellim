<<<<<<< HEAD
function [particleNData] = ParticleData(tracks, k)
close all

particleNData = tracks(tracks(:,end)==k, :);

figure()
plot(particleNData(:,18),movmean(particleNData(:,5),2))

max(particleNData(:,5))
=======
function [particleNData] = ParticleData(tracks, k)
close all

particleNData = tracks(tracks(:,end)==k, :);

figure()
plot(particleNData(:,18),movmean(particleNData(:,5),2))

max(particleNData(:,5))
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
