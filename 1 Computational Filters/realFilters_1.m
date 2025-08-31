% This script:
% 1) Reads in the transmittance spectra of the 16 bandpass filters
% 2) Calculates their effective spectral sensitivities by multiplying them by the spectrally non-neutral components of the camera (IR filter, lens, camera sensor), which are loaded by the script from
%     the file 'cameraComponentsSpectra.mat'. 
% 3) outputs these effective spectral sensitivities in a .mat file called 'realFilters.mat'
% 4) plots these effective spectral sensitivities and saves this as a jpeg file named 'realFilters.jpeg'
% The script should be run from within the '1 Computational Filters/' directory.

clear
close all
load('cameraComponentsSpectra.mat')

files = dir('Transmission*.txt');
numfiles = length(files);
filters = zeros(size(camera,1),numfiles);
colors = rand(numfiles,3);
lambda = (270:991); 
for i = 1:numfiles
    fID = fopen(files(i).name);
    textscan(fID,'%f %f','HeaderLines',100);
    data = textscan(fID,'%f %f',960);
    fclose(fID);
    y = interp1(data{1,1}, data{1,2}, lambda).';
    y = y / 100;
    y = y .* camera.IRblockingFilter .* camera.CoastalOptLens .* camera.JAIsensor;     
    filters(:,i) = y;
    plot(lambda,y, 'LineWidth', 1.5, 'Color', colors(i,:))
    xlim([270 991])
    ylim([0 1])
    hold on
end

ylabel('relative sensitivity')
xlabel('wavelength')
yticks([])
set(gcf, 'Position', [1 1 900 350]) 

print('realFilters.jpeg', '-r600', '-djpeg')
save('realFilters.mat', 'filters')