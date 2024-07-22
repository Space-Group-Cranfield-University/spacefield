% This script provides an example of uniform sampling from a
% two-dimensional ring.

clear
close all

Options.R_1 = 4000;
Options.R_2 = 8000;
Options.N = 10000;
Options.dim = "2D";

rSample = initialiseSampleShell(Options);

figure()
scatter(rSample(:,1), rSample(:,2))
axis equal