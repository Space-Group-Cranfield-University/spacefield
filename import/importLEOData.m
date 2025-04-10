%% Import SpaceTrack catalogue data

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 24);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["INTLDES", "NORAD_CAT_ID", "OBJECT_TYPE", "SATNAME", "COUNTRY", "LAUNCH", "SITE", "DECAY", "PERIOD", "INCLINATION", "APOGEE", "PERIGEE", "COMMENT", "COMMENTCODE", "RCSVALUE", "RCS_SIZE", "FILE", "LAUNCH_YEAR", "LAUNCH_NUM", "LAUNCH_PIECE", "CURRENT", "OBJECT_NAME", "OBJECT_ID", "OBJECT_NUMBER"];
opts.VariableTypes = ["string", "double", "categorical", "categorical", "categorical", "datetime", "categorical", "string", "double", "double", "double", "double", "string", "double", "double", "categorical", "double", "double", "double", "string", "categorical", "categorical", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["INTLDES", "DECAY", "COMMENT", "LAUNCH_PIECE", "OBJECT_ID"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["INTLDES", "OBJECT_TYPE", "SATNAME", "COUNTRY", "SITE", "DECAY", "COMMENT", "RCS_SIZE", "LAUNCH_PIECE", "CURRENT", "OBJECT_NAME", "OBJECT_ID"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "LAUNCH", "InputFormat", "dd/MM/yyyy", "DatetimeFormat", "preserveinput");

% Import the data
SATCAT = readtable(importDirectory+"\satcat-undecayed-leo-no-starlink.csv", opts);


%% Clear temporary variables
clear opts