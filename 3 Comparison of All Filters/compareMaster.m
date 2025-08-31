% For the bird specimen data, this script can optionally be used to run 'compareCustomComputStatNdim.m' and 
% 'plotAllBirdDataNDim.m' within each bird specimen subfolder. If removeDark = 1 here and in 'compareCustomComputStatNdim',
% the script will also output the total number of patches analyzed, the number discarded because they were below the 
% threshold darkness level set in 'compareCustomComputStatNdim', as well as the number discarded expressed as
% a percentage of total patches.
clear
%% user-adjustable parameters
removeDark = 0;
addpath('/Users/baw1743/COCKATOO/16 Filters Data/EOvSpectrocam/Files for Publication/3 Comparison of All Filters')
%%
currentDir = pwd;
files = dir(currentDir);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
for i = 3:numel(subFolders)
    cd([subFolders(i).folder, '/', subFolders(i).name])
    R2 = compareCustomComputStatNdim;
    cd ..
end

R2 = plotAllBirdDataNDim;

if removeDark == 1
    totalPatches = 0;
    totalDiscarded = 0;
    for i = 3:numel(subFolders)
        load([subFolders(i).folder, '/', subFolders(i).name, '/removeDarkInfo.mat'])
        totalPatches = totalPatches + numPatches;
        totalDiscarded = totalDiscarded + numDiscarded;
    end
    
    totalPatches
    totalDiscarded
    percDiscarded = totalDiscarded / totalPatches
    save('removeDarkInfo.mat', 'totalPatches', 'totalDiscarded', 'percDiscarded')
end