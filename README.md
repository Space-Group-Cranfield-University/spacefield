# Spacefield
A MATLAB library for simulating constellations observing space debris, developed at Cranfield University.

# What can you do with the library
Currently, the library contains the necessary code to replicate the results in the 2025 paper "Preliminary analysis and design of an optical space surveillance and tracking constellation for LEO coverage", available [here](https://doi.org/10.1016/j.actaastro.2025.02.019) as an open-access article. The paper discusses how to decouple optical sensor and constellation design, and how space-based sensors can be optimally deployed to maximise coverage of targets over the LEO region. The application domain is that of the observation of space debris, to improve the sustainability of the space enviornment by gathering better data.

All the scripts required to replicate the results of the paper can be found in the <kbd>examples</kbd> folder, where further information is available.

# Installation
1. Download the library as a ZIP file to a local folder.
2. Extract the ZIP file content.
3. Open MATLAB and select the local folder as the main path of the workspace.
4. Add the path to MATLAB, for instance using <kbd>addpath(genpath("spacefield-main"))</kbd> at the beginning of your script.
