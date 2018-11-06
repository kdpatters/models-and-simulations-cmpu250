%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Air_molecule definition
classdef Air_molecule
   properties(Constant = true)
       ambient_ke = 50; % KE_per_A of particles entering simulation
       max_ke = 200; % At what KE should red-hot be?  Only has cosmetic
                     % effect...
       % Alpha can be on [0, pi).  If it's smaller, some particles will
       % respawn on edge travelling parallel to edge
       alpha = 0; % Constant to shrink interval for angle of velocity
                       % air molecules can have when they respawn
   end
   properties
      collisions_enabled = true;
      pos = ones(1, 2, 'double'); % Position
      rpos = ones(1, 2, 'double'); % Rounded position
      vel = zeros(1, 2, 'double'); % Constant velocity
      ke; % Kinetic energy
      T = 50; % Temperature (KE / A)
      handle = -1; % Handle in particles array
      marker = 'o';
      mass = 1; % Mass
      mobile = true; % Is particle movable?
      width = 2;
      height = 2;
      area; % Area of particle; dependent property to `width` and `height`
      therm_cond = 0.5; % Rate of heat transfer; must bet [0, 1]
   end
   methods
      % Calculate kinetic energy per unit area
      function ke = get.ke(obj)
          ke = obj.T * obj.area;
      end
      % Calculate area
      function A = get.area(obj)
        A = abs(obj.width * obj.height);             
      end
      function obj = particle(position, velocity, T)
         % Constructor either takes all initial conditions as input
         % or none, in which case the properties are their defaults.
         % :param position: list with x and y position as floats
         % :param velocity: float for initial velocity
         % :param T: float for kinetic energy of particle per area

         % Parameters are given for particle
         if nargin > 0
             obj.pos = position;
             obj.vel = velocity;
             obj.T = T;
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
          % Randomly sets the velocity of the particle: uses temperature 
          % `T` for speed and picks random direction.
          
          % Assign new velocities
          speed = obj.T;
          random_vector = ((2 * rand(1, length(obj.vel))) - 1);
          unit_vector = random_vector ./ norm(random_vector);
          obj.vel = speed * unit_vector;
      end
      function obj = respawn(obj, grid, x_bounds, y_bounds, curr_step, ...
              n_steps, dt)
          % Resets `Air_molecule` to ambient kinetic energy levels, set
          % position randomly to one of the four edges of the
          % simulation, and give `Air_molecule` a random velocity +/- (pi/2
          % - obj.alpha)
          % of the normal for that edge of the simulation, returning new
          % object.
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Air_molecule objects
          %     in `particles` array
          % :param x_bounds: array containing minimum and maximum x
          %     coordinates between which `Air_molecule` can respawn
          % :param x_bounds: array containing minimum and maximum y
          %     coordinates between which `Air_molecule` can respawn
          % :param curr_step: number of steps elapsed in the direction of
          %     the velocity vector
          % :param n_steps: total number of steps Air_molecule will need to
          %     travel in direction of its velocity vector
          % :param dt: float representing size of timestep being evaluated
          
          obj.T = obj.ambient_ke; % Reset kinetic energy to ambient
                        
          % Choose random position along 4 edges of simulation
          edge = randi([1 4], 1, 1);
          if (edge == 1) % Left edge of canvas
              obj.pos = [x_bounds(1), randi(y_bounds, 1)];
              r_pos = round(obj.pos);
              shift = 3 * pi / 2;
          elseif (edge == 2) % Right edge of canvas
              obj.pos = [x_bounds(2), randi(y_bounds, 1)];
              r_pos = round(obj.pos);
              shift = pi / 2;
          elseif (edge == 3) % Top edge of canvas
              obj.pos = [randi(x_bounds, 1), y_bounds(2)];
              r_pos = round(obj.pos);
              shift = pi;
          else % Bottom edge of canvas
              obj.pos = [randi(x_bounds, 1), y_bounds(1)];
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
          theta = shift + (pi - 2 * obj.alpha) * rand() + obj.alpha;
          unit_vector = [cos(theta) sin(theta)];
          speed = obj.T;
          obj.vel = speed * unit_vector;
          
          % Move Air_molecule the rest of its steps?
          % I've skipped this code for now since implementing a recursive
          % call to the move function could result in the simulation
          % getting stuck in an infinitely loop with a high enough density
          % of particles.
      end
      function [particles, obj, cpu, wall] = ...
              move(obj, particles, dt, grid, cpu, maxes, wall)
          % Move Air_molecule according to time step size `dt`, current
          % position, and current velocity (based on `KE_per_A` property)
          % :param particles: array containing the air molecule 
          %     `Air_molecule`s
          % :param dt: float representing size of timestep being evaluated
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Air_molecule objects
          %     in `particles` array
          % :param cpu: array containing CPU `Heatspreader` objects
          % :param maxes: vector containing max x and y coordinate in grid
          % :param wall: array containing wall `Case_particle` objects
          
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
              % current position before Air_molecule reaches final 
              % destination
              n_steps = dist;
              small_dt = dt / n_steps; % Time between each smaller step
              curr_step = 0;
              path_pos = 0;
              
              % Air_molecule hasn't collided with anything yet
              in_bounds = true;
              no_collisions = true;
              
              while and(no_collisions, ...
                      and(in_bounds, path_pos < final_pos))
                  curr_step = curr_step + 1;
                  path_pos = obj.pos + direction; % Find new pos on path
                  rpath_pos = round(path_pos); % Find rounded position
                  
                  % Try to move one unit using rounded position
                  [obj, no_collisions, in_bounds, particles, cpu, wall] = ...
                      obj.move_helper(rpath_pos, particles, ...
                      small_dt, grid, cpu, curr_step, n_steps, ...
                      maxes, wall);
                  
                   if not(no_collisions) % Particle collided with something
                       % Recursively call move function
                       new_dt = dt - n_steps * small_dt;
                       [particles, obj, cpu] = ...
                        move(obj, particles, new_dt, grid, cpu, ...
                           maxes, wall);
                      
                   elseif and(in_bounds, no_collisions)
                    obj.pos = path_pos; % Update position of `Air_molecule`
                  end
              end
          end
      end
      function [obj, no_collisions, in_bounds, particles, cpu, wall] = ...
              move_helper(obj, new_pos, particles, dt, grid, cpu, ...
              curr_step, n_steps, maxes, wall)
          % Given a new potential position, tries to move object along its
          % velocity vector one unit and returns that object,
          % :param new_pos: integer vector with approximate new position
          % :param particles: array containing the air molecule 
          %     `Air_molecule` objects
          % :param dt: float representing size of timestep being evaluated
          % :param grid: 2D table containing gridcell objects with
          %     properties of `empty` and references to Air_molecule objects
          %     in `particles` array
          % :param cpu: array containing CPU `Heatspreader` objects
          % :param curr_step: number of steps elapsed in the direction of
          %     the velocity vector
          % :param n_steps: total number of steps Air_molecule will need to
          %     travel in direction of its velocity vector
          % :param maxes: vector containing max x and y coordinate in grid
          % :param wall: array containing wall `Case_particle` objects
          
          no_collisions = true; % Set default of no_collisions
          in_bounds = true; % Set default
          
          % Check if position of particle will even change
          if (obj.rpos == new_pos) % Don't do anything; particle won't move
              cpu_h = grid(new_pos(1), new_pos(2)).cpu_h;
              
              if not(cpu_h == -1) % Then exchange energy with CPU
                victum = cpu(cpu_h);
                [obj, victum] = obj.exchange_energies(victum);
                % Update CPUs array
                cpu(cpu_h) = victum;
              end
              
          % Check if particle will go out of the simulated area
          elseif or(or((new_pos(1) <= 0), (new_pos(1) > maxes(1))), ...
                  or((new_pos(2) <= 0), (new_pos(2) > maxes(2))))
              
              % Set bounds toggle
              in_bounds = false;
              
              % Replace particle by another random particle enering
              % simulated area from random direction
              x_bounds = [1 maxes(1)];
              y_bounds = [1 maxes(2)];
              obj = obj.respawn(grid, x_bounds, y_bounds, curr_step, ...
                 n_steps, dt);
          else
              
              % New positon is empty
              if grid(new_pos(1), new_pos(2)).empty
                  obj.rpos = new_pos; % Accept new position
                  
              else
                  wall_h = grid(new_pos(1), new_pos(2)).wall_h;
                  
                  if not(wall_h == -1) % Particle hit a wall
                      
                      [obj, wall(wall_h)] = ...
                          obj.exchange_energies(wall(wall_h));
                      
                       % Particle bounces off wall
                      [obj, particles, cpu, wall] = ...
                          obj.bounce(particles, cpu, wall, wall_h);
                          no_collisions = false;
                  else % Collision!
                      victum = particles(grid(new_pos(1), ...
                          new_pos(2)).particle_h);

                      if and(obj.collisions_enabled, ...
                              victum.collisions_enabled)
                          % Set collisions toggle
                          no_collisions = false;
                          % Recalculate properties based on results 
                          % of collision
                          [obj, victum] = obj.collision(victum);
                          % Update particles array
                          particles(grid(new_pos(1), ...
                              new_pos(2)).particle_h) = ...
                              victum;

                      else
                          obj.rpos = new_pos; % Accept new position
                      end
                  end
              end
              
              % Air_molecule being hit
              cpu_h = grid(new_pos(1), new_pos(2)).cpu_h;
              
              if not(cpu_h == -1) % Then exchange energy with CPU
                  victum = cpu(cpu_h);
                  [obj, victum] = obj.exchange_energies(victum);
                  
                  % Update CPUs array
                  cpu(cpu_h) = victum;
              end
          end
      end
      function [obj, particles, cpu, wall] = ...
            bounce(obj, particles, cpu, wall, wall_h)
          % Give particle random velocity after if bounces a solid 
          % part of the computer case
          % :param cpu: array containing CPU `Heatspreader` objects
          % :param wall: array containing wall `Case_particle` objects
          
          % Wall is horizontal
          if (wall(wall_h).width > wall(wall_h).height)
              % Check direciton molecule was travelling
              if (obj.vel(2) > 0) % Travelling in positive y
                  shift = pi;
              else % Travelling in negative y
                  shift = 0;
                  
              end
          % Wall is vertical
          else
              % Check direciton molecule was travelling
              if (obj.vel(1) > 0) % Travelling in positive x
                  shift = pi / 2;
              else % Travelling in negative x
                  shift = 3 * pi / 2;
              end
          end
          
          % Choose choose random velocity within +/-(pi/2 - alpha) of 
          % the normal to given edge particle will spawn on
          theta = shift + (pi - 2 * obj.alpha) * rand() + obj.alpha;
          unit_vector = [cos(theta) sin(theta)];
          speed = obj.T;
          obj.vel = speed * unit_vector;
      end
      function [obj, p2] = exchange_energies(obj, p2)
          % Exchange amounts of kinetic energy based on size of each and
          % how much time has elapsed, returning changed objects
          % :param obj: a particle object
          % :param p2: a particle object
          % :param dt: float representing size of timestep being evaluated
          
          % Get velocities of particles
          ke1 = obj.ke;
          A1 = obj.area;
          ke2 = p2.ke;
          A2 = p2.area;
          
          % Weight by size of mediums (since size is proportional to mass)
          avg_ke_per_A = (ke1 + ke2) / (A1 + A2);
          
          % Weight change by time difference
          tc1 = obj.therm_cond;
          tc2 = p2.therm_cond;
          tc_avg = sqrt(tc1 * tc2); % Average thermal conductivity
          delta_ke1 = avg_ke_per_A - obj.T;
          delta_ke2 = avg_ke_per_A - p2.T;
          
          % Change in `T` based on thermal conductivity
          obj.T = obj.T + delta_ke1 * tc_avg;
          p2.T = p2.T + delta_ke2 * tc_avg;
      end
      function [obj, p2] = collision(obj, p2)
          % Recalculate velocities after collision between two objects and
          % return both objects.  Formula used is derived from conservation
          % of momentum and conservation of kinetic energy (velocity) in
          % elastic collisions.
          % :param p2: Another `Air_molecule` object
          
          % Get velocities of particles
          v1 = obj.vel;
          v2 = p2.vel;
          m1 = obj.mass;
          m2 = p2.mass;
          
          % Formulae for velocities after elastic collision
          u1 = (v1 * (m1 - m2) + 2 * m2 * v2) / ...
              (m1 + m2);
          u2 = (v2 * (m2 - m1) + 2 * m1 * v1) / ...
              (m1 + m2);
          
          % Return two objects with new properties from after collision
          obj.vel = u1;
          p2.vel = u2;
          obj.T = norm(u1);
          p2.T = norm(u2);
      end
   end
end