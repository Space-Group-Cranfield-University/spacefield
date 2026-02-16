% initialiseSampleShell samples uniformly from a shell, i.e. a volume
% contained between two spheres sharing the same centre. The collection of 
% point can be used for Monte Carlo (MC) integration. The sample points are
% identified by their cartesian coordinates.
%
%   rSample = initialiseSampleShell(Options)
%
%   Inputs:
%   -   Options : structure containing the following fields:
%       .R_1    : lower radius of the shell (0 in the case of a sphere)
%       .R_2    : higher radius of the shell
%       .N      : number of Monte Carlo points
%       .dim    : dimension. If it is "2D", then the sampling is
%               generated over a two-dimensional ring instead of a 3D
%               shell.
%
%   Outputs:
%   -   rSample : generated points in Cartesian coordinates

function rSample = initialiseSampleShell(Options)

    if isfield(Options, "dim") == 0
        Options.dim = "3D";
    end

    switch strcmp(Options.dim, "2D")
        case 1
            rSample = sample2D(Options);
        otherwise
            rSample = sample3D(Options);
    end

end

function rSample = sample3D(Options)

    N = Options.N;
    R_1 = Options.R_1;
    R_2 = Options.R_2;

    dR3 = ( R_2^3 - R_1^3 );
    R13 = R_1^3;
    twoPi = 2 * pi;

    r_sph = zeros(N, 3);
    rSample = zeros(N, 3);

    for k = 1:N
        rand_r = rand();
        rand_t = rand();
        rand_p = rand();

        r_temp      = ( rand_r * dR3 + R13 )^(1/3);
        theta_temp  = 2 * asin( sqrt( rand_t ) );
        psi_temp    = twoPi * rand_p;

        r_sph(k, :) = [ r_temp, theta_temp, psi_temp ];

        rSample(k, :) = [ r_temp * sin( theta_temp ) * cos(psi_temp), ...
                        r_temp * sin( theta_temp ) * sin(psi_temp), ...
                        r_temp * cos( theta_temp ) ]; 
    end

end

function rSample = sample2D(Options)

    N = Options.N;
    R_1 = Options.R_1;
    R_2 = Options.R_2;

    dR2 = ( R_2^2 - R_1^2 );
    R12 = R_1^2;
    twoPi = 2 * pi;

    r_sph = zeros(N, 2);
    rSample = zeros(N, 2);

    for k = 1:N
        rand_r = rand();
        rand_t = rand();

        r_temp      = ( rand_r * dR2 + R12 )^(1/2);
        theta_temp  = twoPi * rand_t;

        r_sph(k, :) = [ r_temp, theta_temp ];

        rSample(k, :) = [ r_temp * cos( theta_temp ), r_temp * sin( theta_temp )]; 
    end

end