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
classdef Particle
   properties(Constant = true)
       ambient_ke = 50;
       max_ke = 200;
   end
   properties
      collisions_enabled = true;
      pos = ones(1, 2, 'double'); % Position
      rpos = ones(1, 2, 'double'); % Rounded position
      vel = zeros(1, 2, 'double'); % Constant velocity
      ke; % Essentially speed; dependent property to `ke_per_A` and `area`
      ke_per_A = 50; % Density of kinetic energy per area
      handle = -1; % Handle in particles array
      marker = 'o';
      mass = 1; % Mass
      mobile = true; % Is particle movable?
      width = 1;
      height = 1;
      area; % Area of particle; dependent property to `width` and `height`
   end
   methods
      % Calculate kinetic energy per unit area
      function ke = get.ke(obj)
          ke = obj.ke_per_A * obj.area;
      end
      % Calculate area
      function A = get.area(obj)
        A = obj.width * obj.height;             
      end
      function obj = particle(position, velocity, ke_per_A)
         % Constructor either takes all initial conditions as input
         % or none, in which case the properties are their defaults.
         % :param position: list with x and y position as floats
         % :param velocity: float for initial velocity
         % :param ke_per_A: float for kinetic energy of particle per area

         % Parameters are given for particle
         if nargin > 0
             obj.pos = position;
             obj.vel = velocity;
             obj.ke_per_A = ke_per_A;
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
          speed = obj.ke_per_A;
          random_vector = ((2 * rand(1, length(obj.vel))) - 1);
          unit_vector = random_vector ./ norm(random_vector);
          obj.vel = speed * unit_vector;
      end
      function obj = respawn(obj, grid, x_bounds, y_bounds, curr_step, ...
              n_steps, dt)
          % Resets `Particle` to ambient kinetic energy levels, set
          % position randomly to one of the four edges of the
          % simulation, and give `Particle` a random velocity +/- pi/4
          % of the normal for that edge of the simulation, returning new
          % object.
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Particle objects
          %     in `particles` array
          % :param x_bounds: array containing minimum and maximum x
          %     coordinates between which `Particle` can respawn
          % :param x_bounds: array containing minimum and maximum y
          %     coordinates between which `Particle` can respawn
          % :param curr_step: number of steps elapsed in the direction of
          %     the velocity vector
          % :param n_steps: total number of steps Particle will need to
          %     travel in direction of its velocity vector
          % :param dt: float representing size of timestep being evaluated
          
          obj.ke_per_A = obj.ambient_ke; % Reset kinetic energy to ambient
                        
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
              if (curr_step < n_steps)
                curr_step = curr_step + 1; % Progress one step
              end
              obj = obj.respawn(grid, x_bounds, y_bounds, curr_step, ...
              n_steps, dt);
          else
              obj.rpos = r_pos;
          end
          
          % Choose choose random velocity within +/-(pi/2 - alpha) of 
          % the normal to given edge particle will spawn on
          alpha = pi / 4; % Constant to shrink interval by so no
                          % particles spawn with velocity parallel to edge
          theta = shift + (pi - 2 * alpha) * rand() + alpha;
          unit_vector = [cos(theta) sin(theta)];
          speed = obj.ke_per_A;
          obj.vel = speed * unit_vector;
          
          % Move Particle the rest of its steps?
          % I've skipped this code for now since implementing a recursive
          % call to the move function could result in the simulation
          % getting stuck in an infinitely loop with a high enough density
          % of particles.
      end
      function [particles, obj, cpu] = ...
              move(obj, particles, dt, grid, cpu)
          % Move Particle according to time step size `dt`, current
          % position, and current velocity (based on `KE_per_A` property)
          % :param particles: array containing the air molecule `Particle`s
          % :param dt: float representing size of timestep being evaluated
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Particle objects
          %     in `particles` array
          % :param cpu: array containing CPU `Heatspreader` objects
          
          % Get new position of particle after one time step
          final_pos = obj.pos + obj.vel * dt;
          % Check rounded position to see if particle is moving into
          % occupied grid cell
          r_pos = round(final_pos); % Rounded position
          
          % Check if position of particle will even change
          if (obj.rpos == r_pos)
              obj.pos = final_pos;
              % Rounded position is the same so we don't need to change it
          else % Calculate path of particle
              % Find distance of travel
              dist = norm(obj.pos - final_pos);
              % Find direction of travel
              direction = obj.vel / norm(obj.vel);
              % Number of iterations of adding direction unit vector to
              % current position before Particle reaches final destination
              n_steps = dist;
              small_dt = dt / n_steps; % Time between each smaller step
              curr_step = 0;
              % Particle hasn't collided with anything yet
              no_collisions = true;
              while and(no_collisions, not(obj.pos == final_pos))
                  curr_step = curr_step + 1;
                  path_pos = obj.pos + direction; % Find new pos on path
                  rpath_pos = round(path_pos); % Find rounded position
                  % Test if rounded position will create a collision
                  [obj, no_collisions, particles, cpu] = ...
                      obj.move_helper(obj, rpath_pos, particles, ...
                      small_dt, grid, cpu, curr_step, n_steps);
                  if no_collisions
                    obj.pos = path_pos; % Update position of `Particle`
                  end
              end
          end
      end
   end
   methods(Static = true)
      function [p, no_collisions, particles, cpu] = ...
              move_helper(p, new_pos, particles, dt, grid, cpu, ...
              curr_step, n_steps)
          % Given a new potential position, tries to move object along its
          % velocity vector one unit and returns that object,
          % :param p: Particle object
          % :param new_pos: integer vector with approximate new position
          % :param particles: array containing the air molecule `Particle`
          %     objects
          % :param dt: float representing size of timestep being evaluated
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Particle objects
          %     in `particles` array
          % :param cpu: array containing CPU `Heatspreader` objects
          % :param curr_step: number of steps elapsed in the direction of
          %     the velocity vector
          % :param n_steps: total number of steps Particle will need to
          %     travel in direction of its velocity vector
          
          no_collisions = true; % Set default of no_collisions
          [max_y, max_x] = size(grid);
          
          in_bounds = true;
          % Check if position of particle will even change
          if (p.rpos == new_pos) % Don't do anything; particle won't move
          % Check if particle will go out of the simulated area
          elseif or(or((new_pos(1) <= 0), (new_pos(1) > max_x)), ...
                  or((new_pos(2) <= 0), (new_pos(2) > max_y)))
              % Set collisions toggle
              no_collisions = false;
              in_bounds = false;
              % Replace particle by another random particle enering
              % simulated area from random direction
              x_bounds = [1 max_x];
              y_bounds = [1 max_y];
              p = p.respawn(grid, x_bounds, y_bounds, curr_step, ...
                  n_steps, dt);
          % New positon is empty
          elseif grid(new_pos(1), new_pos(2)).empty
              p.rpos = new_pos; % Accept new position
          else % Collision!
                victum = particles(grid(new_pos(1), new_pos(2)).particle_h);
                if and(p.collisions_enabled, victum.collisions_enabled)
                    % Set collisions toggle
                    no_collisions = false;
                    % Recalculate properties based on results of collision
                    [victum, p] = p.collision(victum, p, dt);
                    % Update particles array
                    particles(grid(new_pos(1), new_pos(2)).particle_h) = ...
                        victum;
                else
                    p.rpos = new_pos; % Accept new position
                end
          end
          if in_bounds
              % Particle being hit
              cpu_h = grid(new_pos(1), new_pos(2)).cpu_h;
              if not(cpu_h == -1) % Exchange energy with CPU
                victum = cpu(cpu_h);
                [victum, p] = p.exchange_energies(victum, p, dt);
                % Update CPUs array
                cpu(cpu_h) = victum;
              end
          end
      end
      function [p1, p2] = exchange_energies(p1, p2, dt)
          % Exchange amounts of kinetic energy based on size of each and
          % how much time has elapsed, returning changed objects
          % :param p1: a particle object
          % :param p2: a particle object
          % :param dt: float representing size of timestep being evaluated
          
          % Get velocities of particles
          ke1 = p1.ke;
          A1 = p1.area;
          ke2 = p2.ke;
          A2 = p2.area;
          
          % Formulae for velocities after elastic collision
          % Weight by size of mediums (since size is proportional to mass)
          avg_ke_per_A = (ke1 + ke2) / (A1 + A2);
          % Weight change by time difference
          p1.ke_per_A = p1.ke_per_A * (1 - dt) + avg_ke_per_A * dt;
          p2.ke_per_A = p2.ke_per_A * (1 - dt) + avg_ke_per_A * dt;
      end
      function [p1, p2] = collision(p1, p2, dt)
          % Get velocities of particles
          v1 = p1.vel;
          v2 = p2.vel;
          m1 = p1.mass;
          m2 = p2.mass;
          
          % Formulae for velocities after elastic collision
          u1 = (v1 * (m1 - m2) + 2 * m2 * v2) / ...
              (m1 + m2);
          u2 = (v2 * (m2 - m1) + 2 * m1 * v1) / ...
              (m1 + m2);
          
          % Return two objects with new properties from after collision
          p1.vel = u1;
          p2.vel = u2;
          
          [p1, p2] = p1.exchange_energies(p1, p2, dt);
      end
   end
end