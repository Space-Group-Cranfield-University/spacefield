function TRG = propagateTargetPopulation(timeVec, TRG, Forces, CONST, PropagationOptions)
    if nargin < 5
        PropagationOptions = odeset('RelTol',1e-10,'AbsTol',1e-10);
    end
    if nargin < 4
        CONST = initializeAstronomicalConstants;
    end
    if nargin < 3
        Forces.Model(1) = "J2";
        Forces.Model(2) = "atmospheric-drag";
        Forces.B = 0.01;
    end
    
    a = buildAccelerationModel(Forces, CONST);
    dynamicsFunction = @(x) [x(4:6, 1); a(x)];
    parfor k = 1:size(TRG, 2)
        [~, xMat] = ode45(@(t, x) dynamicsFunction(x), timeVec, TRG(k).x0, PropagationOptions);
        TRG(k).xMat = xMat;
    end
end