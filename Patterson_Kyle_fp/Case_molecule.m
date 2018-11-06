%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Case_molecule definition
classdef Case_molecule < Air_molecule
    methods
        function obj = particle(mass, position, ke_per_A, tc)
            % Constructor either takes all initial conditions as input
            % or none, in which case the properties are their defaults.
            % :param mass: float for mass of object
            % :param position: list with x and y position as floats
            % :param ke_per_A: float for kinetic energy of particle per A
            % :param tc: float for rate of thermal condictivity

            obj.mobile = false;
            obj.vel = 0;
            
            % Parameters are given for particle
            if nargin > 0
                obj.pos = position;
                obj.ke_per_A = ke_per_A;
                obj.mass = mass;
                obj.thermal_conductivity = tc;
            else
                obj.thermal_conductivity = 1;
            end
        end
    end
end
    