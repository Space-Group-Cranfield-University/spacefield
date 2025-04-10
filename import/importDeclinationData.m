%% Import MASTER declination distribution data for objects above 10 cm

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 19);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ["\t", " "];

% Specify column names and types
opts.VariableNames = ["Var1", "declination", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "objectDensity", "x14", "x15"];
opts.SelectedVariableNames = ["declination", "objectDensity"];
opts.VariableTypes = ["string", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15"], "EmptyFieldRule", "auto");

% Import the data
DeclinationData = readtable(importDirectory+"\distributionDebrisDeclinationProcessed.txt", opts);


%% Clear temporary variables
clear opts