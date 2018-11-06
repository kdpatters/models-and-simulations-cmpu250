%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Heatspreader definition
classdef Heatspreader < Case_molecule
    properties
        heat_gen = 10; % Rate of heat generation (measured in KE A^-1 t)
    end
    methods
        function obj = Heatspreader()
            obj.collisions_enabled = false; % Disable collisions
        end
        function obj = add_heat(obj, dt)
            obj.T = obj.T + ...
                obj.heat_gen * dt; % Add additional heat
        end
    end
end