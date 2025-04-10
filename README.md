# Spacefield
A MATLAB library for simulating constellations observing space debris, developed at Cranfield University.

# What you can do with the library
Currently, the library contains the necessary code to replicate the results in the 2025 paper "Preliminary analysis and design of an optical space surveillance and tracking constellation for LEO coverage", available [here](https://doi.org/10.1016/j.actaastro.2025.02.019) as an open-access article. The paper discusses how to decouple optical sensor and constellation design, and how space-based sensors can be optimally deployed to maximise coverage of targets over the LEO region. The application domain is that of the observation of space debris, to improve the sustainability of the space enviornment by gathering better data.

All the scripts required to replicate the results of the paper can be found in the <kbd>examples</kbd> folder.

The example folder also contains the scripts necessary to replicate the 2025 conference paper "Performance assessment of space-based optical networks for debris observation and tracking in LEO", available [here](https://www.researchgate.net/publication/390627376_Performance_assessment_of_space-based_optical_networks_for_debris_observation_and_tracking_in_LEO). These scripts may require [ESA MASTER](https://sdup.esoc.esa.int/) debris distribution data or [SpaceTrack](https://www.space-track.org/) catalogue data. Example formatting is available in the <kbd>data_example</kbd> folder.

# Installation
1. Download the [last release](https://github.com/Space-Group-Cranfield-University/spacefield/releases) as a ZIP file to a local folder.
2. Extract the ZIP file content.
3. Open MATLAB and select the local folder as the main path of the workspace.
4. Add the path to MATLAB, for instance using <kbd>addpath(genpath("spacefield-main"))</kbd> at the beginning of your script.
