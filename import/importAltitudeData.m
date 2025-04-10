%% Import MASTER altitude distribution data for objects above 10 cm

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["altitude", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "objectDensity", "Var17", "Var18", "Var19", "Var20"];
opts.SelectedVariableNames = ["altitude", "objectDensity"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "double", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var17", "Var18", "Var19", "Var20"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var17", "Var18", "Var19", "Var20"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "altitude", "TrimNonNumeric", true);
opts = setvaropts(opts, "altitude", "ThousandsSeparator", ",");

% Import the data
AltitudeData = readtable(importDirectory+"\distributionDebrisAltitudeProcessed.txt", opts);


%% Clear temporary variables
clear opts