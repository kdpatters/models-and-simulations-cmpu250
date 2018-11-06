%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Particle definition
classdef Particle_part123
   properties(Constant = true)
       ambient_ke = 100;
       mass = 1;
       collisions_enabled = true;
   end
   properties
      pos = ones(1, 2, 'double'); % Position
      rpos = ones(1, 2, 'double'); % Rounded position
      vel = zeros(1, 2, 'double'); % Constant velocity
      mobile = true; % Is particle movable?
      ke = 100; % Essentially speed
      handle = -1; % Handle in particles array
      marker = 'o';
   end
   methods
      function obj = particle(position, velocity, acceleration, ...
              kinetic_energy)
         % Constructor either takes all initial conditions as input
         % or none, in which case the properties are their defaults.
         % :param position: list with x and y position as floats
         % :param velocity: float for initial velocity
         % :param acceleration: float for acceleration
         % :param kinetic energy: float for kinetic energy of particle
         
         % Parameters are given for particle
         if nargin == 3
             obj.pos = position;
             obj.vel = velocity;
             obj.acc = acceleration;
             obj.ke = kinetic_energy;
         end
      end
      function obj = rand_pos(obj, grid, min_coords, max_coords)
          % Randomly sets the coordinates of the particle based on the
          % input min and max coordinates.
          %
          % :param min_coords: float array for minimum values of
          % coordinates that particle can have
          %
          % :param max_coords: float array for maximum values of
          % coordinates that particle can have
          %disp('Generating position...')

          % Iterate through dimensions and assign coordinates
          for i = 1:length(obj.pos)
              obj.pos(i) = randi([min_coords(i) max_coords(i)], 1, 1);
          end
          r_pos = round(obj.pos); % Set rounded position
          
          % Check if grid position is empty
          if (obj.rpos == r_pos)
              % Do nothing since we already set `obj.pos`
              % Rounded position is the same so we don't need to change it
          elseif not(grid(r_pos(1), r_pos(2)).empty)
              % Warning: If grid is densely filled, this recursive call may
              % result in significant lag
              obj = obj.rand_pos(grid, min_coords, max_coords);
          else
              obj.rpos = r_pos;
          end
      end
      function obj = rand_vel(obj)
          % Randomly sets the velocity of the particle: uses kinetic 
          % energy `ke` for speed and picks random direction.
          
          % Assign new velocities
          speed = obj.ke;
          random_vector = ((2 * rand(1, length(obj.vel))) - 1);
          unit_vector = random_vector ./ norm(random_vector);
          obj.vel = speed * unit_vector;
      end
      function obj = respawn(obj, grid, x_bounds, y_bounds)
          %disp('Respawning...')
          
          obj.ke = obj.ambient_ke; % Reset kinetic energy to ambient
                        
          % Choose random position along 4 edges of simulation
          edge = randi([1 4], 1, 1);
          if (edge == 1)
              obj.pos = [x_bounds(1), randi(y_bounds, 1)];
              r_pos = round(obj.pos);
              shift = 3 * pi / 2;
          elseif (edge == 2)
              obj.pos = [x_bounds(2), randi(y_bounds, 1)];
              r_pos = round(obj.pos);
              shift = pi / 2;
          elseif (edge == 3)
              obj.pos = [randi(x_bounds, 1), y_bounds(2)];
              r_pos = round(obj.pos);
              shift = pi;
          else
              obj.pos = [randi(x_bounds, 1), x_bounds(1)];
              r_pos = round(obj.pos);
              shift = 0;
          end
          
          if (obj.rpos == r_pos)
              % Do nothing since we already set `obj.pos`
              % Rounded position is the same so we don't need to change it
          % Check if particle will go out of the simulated area
          % Check if grid position is empty
          elseif not(grid(r_pos(1), r_pos(2)).empty)
              % Warning: If grid is densely filled, this recursive call may
              % result in significant lag
              obj = obj.respawn(grid, x_bounds, y_bounds);
          else
              obj.rpos = r_pos;
          end
          
          % Choose choose random velocity within +/-(pi/2 - alpha) of 
          % the normal to given edge particle will spawn on
          alpha = pi / 4; % Constant to shrink interval by so no
                          % particles spawn with velocity parallel to edge
          theta = shift + (pi - 2 * alpha) * rand() + alpha;
          unit_vector = [cos(theta) sin(theta)];
          speed = obj.ke;
          obj.vel = speed * unit_vector;
      end
      function [particles, obj, grid] = move(obj, particles, dt, grid)
          % :param t: float value for time
          % :param grid: 2D table containing gridcell objects with
          % properties of `empty` and reference to a particle `Particle`
          %disp('Moving...')
          
          % Get new position of particle after one time step
          new_pos = obj.pos + obj.vel * dt;
          % Check rounded position to see if particle is moving into
          % occupied grid cell
          r_pos = round(new_pos); % Rounded position
          [max_y, max_x] = size(grid);
          
          % Check if position of particle will even change
          if (obj.rpos == r_pos)
              %disp('Same pos')
              obj.pos = new_pos;
              % Rounded position is the same so we don't need to change it
          % Check if particle will go out of the simulated area
          % If so, make particle bounce of wall
          elseif (r_pos(1) <= 0)
              obj.vel = [-obj.vel(2) obj.vel(1)];
          elseif (r_pos(1) > max_x)
              obj.vel = [-obj.vel(2) obj.vel(1)];
          elseif (r_pos(2) <= 0)
              obj.vel = [obj.vel(2) -obj.vel(1)];
          elseif (r_pos(2) > max_y)
              obj.vel = [obj.vel(2) -obj.vel(1)];
          elseif grid(r_pos(1), r_pos(2)).empty
              %disp('Empty')
              obj.pos = new_pos; % Accept new position
              obj.rpos = r_pos;
          else % Collision!
             %disp('Collision')
             % Particle being hit
             victum = particles(grid(r_pos(1), r_pos(2)).particle_h);
             
             if and(obj.collisions_enabled, victum.collisions_enabled)
                 % Recalculate properties based on results of collision
                 [victum, obj] = obj.collision(victum, obj);
                 % Update particles array
                 particles(grid(r_pos(1), r_pos(2)).particle_h) = victum;
             end
          end    
      end
   end
   methods(Static = true)
      function [p1, p2] = collision(p1, p2)
          % disp(strcat('Collision at ', int2str(p1.rpos(1)), ...
          %    int2str(p1.rpos(2))));
          
          % Get velocities of particles
          v1 = p1.vel;
          v2 = p2.vel;
          m1 = p1.mass;
          m2 = p2.mass;
          ke1 = p1.ke;
          ke2 = p2.ke;
          
          % Formulae for velocities after elastic collision
          u1 = (v1 * (m1 - m2) + 2 * m2 * v2) / ...
              (m1 + m2);
          u2 = (v2 * (m2 - m1) + 2 * m1 * v1) / ...
              (m1 + m2);
          avg_ke = mean([ke1 ke2]);
          p1.ke = avg_ke;
          p2.ke = avg_ke;
          
          % Return two objects with new properties from after collision
          p1.vel = u1;
          p2.vel = u2;
      end
   end
end