This MATLAB code is part of the thesis on assessing the effects of solar weather radiation effects as constellations proliferate.

The general flow is:
1. StormIsolator_Final.m
2. GeomagneticShielding_Final.m
3. DoseCalculations_Final.m
4. BackgroundRadiation_Final.m
5. (Optional) [SatName]_STK.m
6. StormSeedGenerator_Final.m
7. MonteCarloSimulation_[SatName]_Final.m

Information about neccessary input files, the flow of each individual file, and the organization of files are in each individual code file.

All references are directed to D:\Materials or D:\ so if you download the files to a different directory, these must be updated in each file or the new file path to the download location must be added
--------------

Sample inputs and outputs for each stage are also provided for Flock, SkySat, and WorldView.

For example, the unique result of propagating SPEs to SkySat's orbit is given in "Propagated_Storms_SkySat.mat".

The final dose calculation files before the Monte Carlo portion are named "Storm_Doses_[SatName].mat"

These files are ~100MB each and are therefore hosted in a Google drive at: https://drive.google.com/drive/folders/1GszLt3o7P6S8ve6ze7qPqH6S7giSFBR0?usp=sharing

These files should be downloaded and moved to the "Materials" folder if seeking to replicate results.
--------------

Besides MATLAB, other pre-requisite required programs are:

STK (Coverage) - for calculating revisit rates
Libirbem - for SHIELDOSE-2 calculations
AE9AP9 - for calculating background radiation

--------------

Feel free to post questions or issues and I will do my best to respond promptly!