# photographic-spectrometric-comparison
data & code associated with Tedore C (2024). A comparison of photographic and spectrometric methods to quantify the colours seen by animal eyes. Methods in Ecology and Evolution 15(1), 4–23. https://doi.org/10.1111/2041-210X.14255

Workspace files with pre-aligned images are available at https://doi.org/10.6084/m9.figshare.22284688.v1 (too large for github). 

To re-create the analyses done in the paper, the user must open MATLAB and:
1.  Set ‘1 Computational Filters/’ as the working directory and run the ‘realFilters_1.m’ script, followed by the ‘computFilters_2.m’ script, following any instructions written at the top of the script.
2.  Set ‘2 Statistical Filters/’ as the working directory and:
    a.  Run the ‘StatisticalFilters_photorecMapping_naturalSpectra.m’ script, following the instructions written at the top of the script.
    b.  Run the ‘StatisticalFilters_photorecMapping_colorChart.m’ script, following the instructions written at the top of the script.
3.  Set the working directory to the one containing the images you want to use to compare methods (e.g., 'Pastels Images/' or 'Bird Specimen Images/Alcedo atthis insipida adult male/')
    a.  Run the ‘compareCustomComputStatNdim.m’ function found in the ‘3 Comparison of All Filters’ folder, following the instructions written at the top of the function.
    b.  Repeat (a) on all bird specimen images.
    c.  Set the working directory to the one containing all bird specimen images (i.e., ‘Bird Specimen Images/’). Then run the ‘plotAllBirdDataNDim.m’ function found in the ‘3 Comparison of All Filters’ folder, following the instructions written at the top of the function.
    Note: To complete 3a-c in one step for all species of bird specimens, one can alternatively run the ‘compareMaster.m’ script. Note, however, that there are a number of adjustable parameters that must first be adjusted at the top of the ‘compareCustomComputStatNdim.m’ and ‘plotAllBirdDataNDim.m’ functions.
4.  Set ‘4 MATLAB vs MICA/’ as the working directory and run the ‘MATvMICA_VSbirds.m’ script.

Note that the ‘3 Comparison of All Filters/Pastels Images/’ directory and each bird image directory has the images ordered thusly, which each bandpass filter image in its own folder, and each exposure set of six custom filters in its own folder: 
Bandpass filter 300 (short exposure)
Bandpass filter 300 (long exposure)
Bandpass filter 325 (short exposure)
Bandpass filter 325 (long exposure)
Bandpass filter 350 (short exposure)
Bandpass filter 350 (long exposure)
Bandpass filter 375 (short exposure)
Bandpass filter 375 (long exposure)
Bandpass filter 400 (short exposure)
Bandpass filter 400 (long exposure)
Bandpass filter 425 (short exposure)
Bandpass filter 425 (long exposure)
Bandpass filter 450 (short exposure)
Bandpass filter 450 (long exposure)
Bandpass filter 475 (short exposure)
Bandpass filter 475 (long exposure)
Bandpass filter 500 (short exposure)
Bandpass filter 500 (long exposure)
Bandpass filter 525 (short exposure)
Bandpass filter 525 (long exposure)
Bandpass filter 550 (short exposure)
Bandpass filter 550 (long exposure)
Bandpass filter 575 (short exposure)
Bandpass filter 575 (long exposure)
Bandpass filter 600 (short exposure)
Bandpass filter 600 (long exposure)
Bandpass filter 625 (short exposure)
Bandpass filter 625 (long exposure)
Bandpass filter 650 (short exposure)
Bandpass filter 650 (long exposure)
Bandpass filter 675 (short exposure)
Bandpass filter 675 (long exposure)
Custom-fabricated UVS filter (short exposure)
Custom-fabricated VS filter (short exposure)
Custom-fabricated SWS2 (UVS visual system) filter (short exposure)
Custom-fabricated SWS2 (VS visual system) filter (short exposure)
Custom-fabricated MWS filter (short exposure)
Custom-fabricated LWS filter (short exposure)
Custom-fabricated UVS filter (long exposure)
Custom-fabricated VS filter (long exposure)
Custom-fabricated SWS2 (UVS visual system) filter (long exposure)
Custom-fabricated SWS2 (VS visual system) filter (long exposure)
Custom-fabricated MWS filter (long exposure)
Custom-fabricated LWS filter (long exposure)

The images in the ‘2 Statistical Filters/Pastels Images/’ (duplicates of those in ‘3 Comparison of All Filters/Pastels Images/) and ‘2 Statistical Filters/Spyder Checkr 24 Images/’ directories are ordered the same, but lack the custom avian filter images.
