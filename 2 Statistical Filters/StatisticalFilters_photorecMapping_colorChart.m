%% STATISTICAL FILTERS CONE MAPPING USING A COLOR CHART
% This script:
% 1) reads in the file of the irradiance spectrum to be simulated (Seric xenon lamp)
% 2) reads in the 16 bandpass filter images of the color chart to be used for mapping as well as their metadata files
% 3) subtracts dark noise from the images
% 3) combines the two exposure levels of images taken into HDR images
% 4) aligns the images
% 5) removes under- and over-exposed pixels and replaces them with NaN
% 6) loads the computational filter coefficients file generated for the Nikon camera D7000 camera converted for full-spectrum photography (computed and saved to the appropriate directory by the 'computFilters_2.m' script)
% 7) uses the above file to convert the 16 bandpass filter images into computational Nikon filter images
% 8) has the user interactively select the gray standard and each color patch on the chart, and converts the Nikon filter images to relative quantum catch images by dividing them by the gray standard
% 9) extracts the relative quantum catch of each color patch on the chart 
%10) reads in the reflectance spectra of the color chart
% 11) calculates animal relative quantum catches in response to each color patch on the chart using the above reflectance spectra 
% 12) fits a linear model relating each photoreceptor relative quantum catch to all camera channel relative quantum catches
% 13) saves all linear models in a single .mat file
% Steps 2-4 are skipped if the 'loadWorkspace' parameter is set to 1. This loads a workspace containing previously-aligned images, saving time and energy if the script is to be run multiple times.
% The user must input the desired parameters in the 'User-adjustable parameters' section below (i.e., 'topDir', 'loadWorkspace', 'pastels', 'spectrumUsed', 'animal', 'loadGS', and 'loadSelect').
% The script should be run from within the '2 Statistical Filters/' directory.
% If trainingSet = 'expanded', the script adds simulated chart spectra and camera relative q of different brightnesses to the training set.
clear
%% User-adjustable parameters:
topDir = ''; % specify full path of directory above the '2 Statistical Filters' directory
pastels = 0; % 1 if chart is set of 48 pastels; 0 if chart is 24-color camera calibration chart
loadWorkspace = 1; % when set to 1, loads a workspace containing previously-aligned images ('workspace.mat')
spectrumUsed = 'VIS'; % can be VIS or UVVIS
animal = 'Lcatta'; % use the 'animal' parameter from the 'computFilters_2.m' script
loadGS = 1; % if you already selected the gray standard when running this script previously, it has been saved to file ('grayStandardPosition.mat'). To load this file rather than select it again, set loadGS = 1
loadSelect = 1; % if you already selected the chart color patches when running this script previously, it has been saved to file ('select.mat'). To load this file rather than select them again, set loadSelect = 1
trainingSet = 'expanded'; % can be expanded or orig
factors = [0.001, 0.01, 0.1:0.1:1, 2:12]; % factors by which to multiply copies of training spectra if trainingSet = 'expanded'

%%
animalSensFilename = [topDir, animal, '_idealizedSens300-780.mat']; % idealized animal spectral sensitivities to map to

if pastels == 1
    chartImgsDirect = 'Pastels Images'; % directory where pastel images can be found
    chartReflDirect = [topDir, 'Reflectance Spectra/Pastels Reflectance Spectra Seric']; % directory where reflectance spectra of chart can be found
    outputFilename = [animal, spectrumUsed, 'Mapping_pastels_Seric_', trainingSet, '.mat']; % name of output mapping file
else
    chartImgsDirect = 'Spyder Checkr 24 Images'; % directory where chart images can be found
    chartReflDirect =  [topDir, 'Reflectance Spectra/Spyder Checkr 24 Reflectance Spectra Seric']; % directory where reflectance spectra of chart can be found
    outputFilename = [animal, spectrumUsed, 'Mapping_chart_Seric_', trainingSet, '.mat']; % name of output mapping file
end

if loadWorkspace == 1
    load([chartImgsDirect, '/workspace.mat'])
end

% Computational Filter Coefficient File (for Nikon D7000 with CoastalOpt 60 mm lens):
camCoeffFilename = [topDir, '1 Computational Filters/NikonD7000CoastalOpt60ComputFilterCoeffs.mat'];
grayStandRefl = 0.2; % reflectance of gray standard 
EOfilterPos = 8; % position of blank spot in filterwheel camera allowing photography through bandpass filter mounted on front of lens
numFilters = 16; % number of filters used to build computational Nikon filters
lampSpecFilename = [topDir, 'Seric_AbsoluteIrradiance_12-38-16-170.txt']; % modeled light source 
lambda = (300:700).';
hdrLevels = 2;
    
if loadWorkspace == 0    
    %% read in exposure times
    metaFiles = dir([chartImgsDirect, '/*/Metadata*.txt']);
    exposures = zeros(1, numFilters);
    for f = 1:numFilters*hdrLevels
        fID = fopen([metaFiles(f).folder, '/', metaFiles(f).name]);
        textscan(fID,'%f %f','HeaderLines',17+4*(EOfilterPos - 1));
        metadata = textscan(fID,'%*s %*f %*s %*s %f4.2',1);
        exposures(f) = metadata{1};
        fclose(fID);
    end
    % calculate hdr factors
    hdrLevels = 2; % only works with 2 levels
    for i = 1:hdrLevels
        EOexpos(i,:) = exposures(i:hdrLevels:numFilters*hdrLevels);
    end
    EOhdrFactors = EOexpos(2,:) ./ EOexpos(1,:);
    
    %% read in photos of color chart
    % delete first photograph in every subfolder (bug in camera software results in this image always being overexposed)
    photosToDelete = dir([chartImgsDirect, '/*/Image_I0000*.tiff']);
    for i = 1:numel(photosToDelete)
        delete([photosToDelete(i).folder, '/', photosToDelete(i).name])
    end
    chartPhotos = dir([chartImgsDirect, '/*/*.tiff']);
    
    % preallocation
    darkNoise = zeros(numFilters*hdrLevels,1);
    maxVal = zeros(numFilters*hdrLevels,1);
    minVal = zeros(numFilters*hdrLevels,1);
    numOverExp = zeros(numFilters*hdrLevels,1);
    numUnderExp = zeros(numFilters*hdrLevels,1);
    maxPoss = zeros(numFilters*hdrLevels,1);
    mydata=cell(1,numFilters*hdrLevels);
    
    % read in images, subtract dark noise
    for k = 1:numel(chartPhotos)
        mydata{k} = imread([chartPhotos(k).folder, '/', chartPhotos(k).name]);
        mydata{k} = double(mydata{k})-32768; % converts 16 bit image into 10 bit image
    
        noise = mydata{k}(:,1393:1408);
        darkNoise(k) = ceil(mean2(noise)); % calculates dark noise in each channel 
        mydata{k} = mydata{k}(5:1040,1:1392); % eliminates unused pixels from further processing
        yDim = size(mydata{k},1);
        xDim = size(mydata{k},2);
        numOverExp(k) = sum(mydata{k}(:) == 1023);
        numUnderExp(k) = sum(mydata{k}(:) <= darkNoise(k));
        
        mydata{k} = mydata{k} - darkNoise(k); % subtracts dark noise
        mydata{k}(mydata{k}(:) < 0) = 0; % converts neg numbers (noise) to zero
        
        maxVal(k) = max(max(mydata{k}));
        minVal(k) = min(min(mydata{k}));
        
        maxPoss(k) = 1023 - darkNoise(k);
    end
    
    img = (1:1:numFilters*hdrLevels).';
    diagnosticsTbl = table(img, darkNoise, numUnderExp, numOverExp, minVal, maxVal, maxPoss)
    
    %% CONVERT TO HDR
    if hdrLevels > 1
        % EO filters
        f = 1; % hdr factor index
        for k = hdrLevels:hdrLevels:numFilters*hdrLevels
            if k == 2 && EOhdrFactors(1) == 1 % if was hard to get enough light for BP 300 filter & thus both shots through this filter had the max exposure
                mydata{2} = (mydata{2} + mydata{1}) / 2; % then add these two shots together and divide by 2 (will decrease number of underexposed pixels)
            end
            for i = 1:yDim
                for j = 1:xDim
                    if mydata{k}(i,j) == maxPoss(k) && mydata{k-1}(i,j) < maxPoss(k-1)
                        mydata{k}(i,j) = EOhdrFactors(f) * mydata{k-1}(i,j);
                    end
                end
            end
            f = f + 1;
        end
    end
    mydata = mydata(hdrLevels:hdrLevels:numFilters*hdrLevels);
    maxPoss = maxPoss(hdrLevels:hdrLevels:numFilters*hdrLevels);

    %% WRITE IMAGES TO CHECK FOR OVER- OR UNDER-EXP0SURE IN ROIs: writes grayscale images of each filter image with over- and under-exposed pixels highlighted in purple.
    for k = 1:numFilters
        mydataCheck = cat(3, mydata{k}, mydata{k}, mydata{k});
        mydataCheck = mydataCheck / max(max(max(mydataCheck)));
        for x = 1:yDim
            for y = 1:xDim
                if mydata{k}(x,y) == maxPoss(k) || mydata{k}(x,y) == 0
                    mydataCheck(x,y,:) = [1 0 1];
                end
            end
        end
        filename = sprintf([topDir, '2 Statistical Filters/', chartImgsDirect, '/exposCheck_%d.png'], k);
        imwrite(mydataCheck, filename, 'png');
    end
    
    %% ALIGN IMAGES
    % USING FUNCTION GENERATED BY 'registrationEstimator' app
    % registrationEstimator(mydata{2}, mydata{1}); 
    for i = 1:numFilters-1
        % order: moving, fixed
        if pastels == 1
            mydata{i+1} = registerImages(mydata{i+1}, mydata{i});
        else
            mydata{i+1} = registerImagesSpyder(mydata{i+1}, mydata{i});
        end
        % compareImgs(mydata{i+1}, mydata{i})
    end

    %% REMOVE UNDER- AND OVER-EXPOSED PIXELS 
    % inserts NaN into under- and over-exposed pixels (across all images)
    % For isolated pixels, only works on 1st fixed image from alignment proces. Works for large groups of pixels in other images. 
    % Be sure to visually check for over- or under-exposure in all filter images.
    for k = 1:numFilters
        for x = 1:yDim
            for y = 1:xDim
                if abs(mydata{k}(x,y) - maxPoss(k)) < 0.00001
                    for j = 1:numFilters
                        mydata{j}(x,y) = NaN;
                    end
                elseif abs(mydata{k}(x,y) - 0) < 0.00001 
                    for j = 1:numFilters
                        mydata{j}(x,y) = NaN;
                    end
                end
            end
        end
    end

    %% CHECK THAT IMAGES WELL-ALIGNED
    for k = 1:numFilters
        mydataCheck = mydata{k} / max(max(mydata{k}));
        filename = sprintf([topDir, '2 Statistical Filters/', chartImgsDirect, '/alignCheck_%d.png'], k);
        imwrite(mydataCheck, filename, 'png')
    end

    save([chartImgsDirect, '/workspace.mat'], 'mydata', 'yDim', 'xDim', 'EOexpos')
end

%% CONVERT TO VIRTUAL FILTERS
% load RGB camera virtual filter coefficient file
load(camCoeffFilename); 
% Adjust photos to all have same effective exposures
mydata_exposCorr = cell(1,numFilters);
exposCorr = 2000 ./ EOexpos(2,:);
for i = 1:numFilters
    mydata_exposCorr{i} = mydata{i} * exposCorr(i);
end
% Build virtual filters
numVirtFilts = size(coeffSet,2);
myVirtualData = cell(1,numVirtFilts);
for i = 1:numVirtFilts
    myVirtualData{i} = zeros(yDim, xDim);
    for j = 1:numFilters
        myVirtualData{i} = myVirtualData{i} + coeffSet(j,i) * mydata_exposCorr{j};
    end
end

%% SELECT GRAY STANDARD
if loadGS == 1
    load([chartImgsDirect, '/grayStandardPosition.mat'])
else
    h = figure;
    imgToDispl = cat(3, myVirtualData{1}, myVirtualData{2}, myVirtualData{3});
    for i = 1:3
        chan_i = imgToDispl(:,:,i);
        chan_i = chan_i / mean(chan_i(~isnan(chan_i)));
        imgToDispl(:,:,i) = chan_i ./ (chan_i + 1);
    end
    imgToDispl = imgToDispl / max(max(max(imgToDispl)));
    imshow(imgToDispl);
    message = sprintf('Left click vertices defining gray standard.\nThen double click in the middle to accept it.');
    uiwait(msgbox(message));
    [gsSelect, xi, yi] = roipoly();
    hold on
    plot(xi, yi, 'r-', 'LineWidth', 2);
    pause(0.2)
    close(h)
    save([chartImgsDirect, '/grayStandardPosition.mat'], 'gsSelect', '-mat')   
end

%% WHITE BALANCE
for k = 1:numVirtFilts
    wStand = myVirtualData{k}(gsSelect);
    wStand = wStand(~isnan(wStand));
    wStand = mean(wStand);
    myVirtualData{k} = myVirtualData{k} / wStand;
end

%% SELECT COLOR PATCHES
if loadSelect == 1
    load([chartImgsDirect, '/select.mat'])
    if pastels == 1
        select = cat(3, select(:,:,1:39), select(:,:,43:48)); % removes fluorescent pastels
    end
else
    iter = 1;
    h = figure;
    imgToDispl = cat(3, myVirtualData{1}, myVirtualData{2}, myVirtualData{3});
    imgToDispl = imgToDispl ./ (imgToDispl + 1);
    imshow(imgToDispl);
    message = sprintf('Zoom and pan to see region of interest clearly and then press any key. Single click vertices defining the region of interest.\nThen double click in the middle to accept it. ');
    uiwait(msgbox(message));
    again = true;
    pause
    while again
        [select(:,:,iter), xi, yi] = roipoly();
        iter = iter + 1;
        hold on
        plot(xi, yi, 'r-', 'LineWidth', 2);
        pause
        promptMessage = sprintf('Draw region #%d in the image or Quit?', iter);
        titleBarCaption = 'Continue?';
        button = questdlg(promptMessage, titleBarCaption, 'Draw', 'Quit', 'Draw');
        if strcmpi(button, 'Quit')
            break;
        end
    end
    print([chartImgsDirect, '/select.png'], '-dpng', '-r600')
    close(h)
    save([chartImgsDirect, '/select.mat'], 'select', '-mat')
end

%% CALCULATE CAMERA q
numPatches = size(select, 3);
camera_q = zeros(5,numPatches);
for s = 1:numPatches
    for i = 1:5
        selectedArea = myVirtualData{i}(select(:,:,s));
        camera_q(i,s) = mean(selectedArea(~isnan(selectedArea)));
    end
end
if strcmp(trainingSet, 'expanded')
    camera_q_orig = camera_q;
    camera_q = [];
    for f = 1:length(factors)
        camera_q = cat(2, camera_q, factors(f)*camera_q_orig);
    end
end

%% read in color chart spectra
chartSpectra = zeros(length(lambda), numPatches);
chartFiles = dir([chartReflDirect, '/Reflection*.txt']);
for i = 1:numel(chartFiles)
    fID = fopen([chartFiles(i).folder, '/', chartFiles(i).name]);
    textscan(fID,'%f %f','HeaderLines',100);
    data = textscan(fID,'%f %f',960);
    fclose(fID);
    chartSpectra(:,i) = interp1(data{1,1}, data{1,2}, lambda) / 100;
end
if pastels == 1
    chartSpectra = horzcat(chartSpectra(:,1:39), chartSpectra(:,43:48)); % removes fluorescent pastels
end
if strcmp(trainingSet, 'expanded')
    chartSpectra_orig = chartSpectra;
    chartSpectra = [];
    for f = 1:length(factors)
        chartSpectra = cat(2, chartSpectra, factors(f)*chartSpectra_orig);
    end
    numPatches = size(camera_q, 2);
end

%% read in modeled light source spectrum
fID = fopen(lampSpecFilename); % Seric xenon lamp 
textscan(fID,'%f %f','HeaderLines',100);
data = textscan(fID,'%f %f',960);
fclose(fID);
lambda = (300:700).';
lightSource = interp1(data{1,1}, data{1,2}, lambda);
lightSource = 5.05E15 * lightSource .* lambda * 0.01;

%% calculate animal q
load(animalSensFilename); % animal spectral sensitivities
idealizedSens = idealizedSens(1:401,:);
if strcmp(spectrumUsed, 'VIS') 
    idealizedSens = idealizedSens(101:end,:);
    lightSource = lightSource(101:end);
    chartSpectra = chartSpectra(101:end, :);
end
numPhotorec = size(idealizedSens, 2);
animal_q = zeros(numPhotorec, numPatches);
for i = 1:numPatches
    for p = 1:numPhotorec
        animal_q(p,i) = sum(chartSpectra(:,i) .* lightSource .* idealizedSens(:,p)) / sum(grayStandRefl * lightSource .* idealizedSens(:,p)); % calculates animal relative quantum catch from all training spectra 
    end
end

%% CONE MAPPING
modelData = horzcat(animal_q.', camera_q.');
if numPhotorec == 4
    varNames = {'P1', 'P2', 'P3', 'P4', 'R', 'G', 'B', 'UB', 'UR'};
elseif numPhotorec == 3
    varNames = {'P1', 'P2', 'P3', 'R', 'G', 'B', 'UB', 'UR'};
elseif numPhotorec == 2
    varNames = {'P1', 'P2', 'R', 'G', 'B', 'UB', 'UR'};
end
modelData = array2table(modelData, "VariableNames", varNames); 

% 2-way interactions
if strcmp(spectrumUsed, 'VIS') 
    effects = 'R*G + R*B + G*B';
else
    effects = 'R*G + R*B + R*UB + R*UR + G*B + G*UB + G*UR + B*UB + B*UR + UB*UR';
end
% effects = 'R + G + B + UB + UR'; % without interactions 
lm{1} = fitlm(modelData, ['P1 ~ ', effects]);
lm{2} = fitlm(modelData, ['P2 ~ ', effects]);
if numPhotorec > 2
    lm{3} = fitlm(modelData, ['P3 ~ ', effects]);
end
if numPhotorec > 3
    lm{4} = fitlm(modelData, ['P4 ~ ', effects]);
end
save(outputFilename, 'lm')