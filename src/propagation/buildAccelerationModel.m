function a = buildAccelerationModel(Forces, CONST)
    if nargin < 2
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 1
        Forces.Model = ["J2", "atmospheric-drag"];
        Forces.B = 0.01;
        Forces.densityModel = "NRLMSIS-00";
    end
    if ~isfield(Forces, "densityModel")
        Forces.densityModel = "NRLMSIS-00";
    end
    a = @(x) getKeplerianAcceleration(x, CONST);
    if isfield(Forces, "Model")
        for k = 1:size(Forces.Model, 2)
            if strcmp(Forces.Model(k), "J2")
                a = @(x) a(x) + getPerturbationJ2(x, CONST);
            end
            if strcmp(Forces.Model(k), "atmospheric-drag")
                a = @(x) a(x) + getPerturbationDrag(x, Forces.B, Forces.densityModel, CONST);
            end
        end
    end
end