%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Grid cell definition
classdef Gridcell
    properties
      empty = true; % Is cell empty?
      particle_h = -1; % Reference to air molecule in grid cell
      wall_h = -1; % Reference to case molecule in grid cell
      cpu_h = -1; % Reference to cpu in grid cell
   end
   methods
       function obj = Gridcell(my_particle)
           % :param my_particle: particle object
           if (nargin > 0)
               obj.particle_h = my_particle;
               obj.empty = false;
           end
       end
   end
end
           