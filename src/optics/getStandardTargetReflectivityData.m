function TRG_DATA = getStandardTargetReflectivityData()
    % getStandardTargetReflectivityData() gets standard space object 
    % reflectivity data.
    %
    % SYNTAX:
    %   TRG_DATA = getStandardTargetReflectivityData()
    %
    % INPUTS:
    %   - None
    % OUTPUTS:
    %   - TRG_DATA  : Data structure containing target reflectivity data
    %       > TRG_DATA.rho_B    : Bond albedo
    %       > TRG_DATA.k_m      : Diffuse / specular mixing coefficient
    %
    % Gets standard space object reflectivity data in a data structure with 
    % fields:
    %   .rho_B  : Bond albedo, assumed to be 0.2
    %   .k_m    : Diffuse / specular mixing coefficient, assumed to be 0.25
    %
    % Note: a low k_m means an almost perfect specular reflection
    %
    % References:
    %   - Clemens SD. On-Orbit Resident Space Object (RSO) Detection Using
    %   Commercial Grade Star Trackers; 2019
    %   - Spiller D, Magionami E, Schiattarella V, Curti F, Facchinetti C, 
    %   Ansalone L, et al. On-orbit recognition of resident space objects 
    %   by using star trackers. Acta Astronautica. 2020 Dec;177:478-96
    %
    % Written by: Antonio D'Anniballe
    % © Cranfield University, 2026

    TRG_DATA.rho_B = 0.2; % Target reflectivity (Bond albedo)
    TRG_DATA.k_m = 0.25; % Diffuse / specular mixing coefficient
end