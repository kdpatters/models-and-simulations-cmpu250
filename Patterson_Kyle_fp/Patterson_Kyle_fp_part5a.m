%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (5)
clear;
part_name = 'part5';

%%% Parameters

% Output parameters
% Size of simulation environment (points/pixels)
canvas_width = 100;
half_width = canvas_width / 2;
canvas_height = 100;
half_height = canvas_height / 2;

report_interval = 0.1; % Spacing between text output reports (seconds)
img_interval = 0.01; % Spacing between image outputs for GIF (seconds)

end_time = 60.0; % Length of simulation (seconds)
dt = 0.01; % Timestep (seconds)
num_points = ceil(end_time / dt) + 1;
% How fast time should change in the simulation relative to real time
playback_speed = 1;
pause_time = dt / playback_speed;

% Ambient temperature
temp = Air_molecule;
ambient_ke = temp.ambient_ke;
max_ke = temp.max_ke;

% Initial numbers of particles
% Number of particles inside case
init_n_inside = 30;

% Number of particles outside of case
init_n_outside = 30;

% Constant increase of overall KE on CPU per second

% Test different `ke_inc`
hold on;
title('CPU Temperature Over Time')
xlabel('Time (seconds)');
ylabel('Kinetic energy');
ke_inc_list = 2:3:8;
for ke_inc = ke_inc_list

    total_particles = init_n_inside + init_n_outside;
    particles(total_particles) = Air_molecule; % Initialize empty object 
                                               % array

    % Create table to store positions of all particles
    grid(canvas_width, canvas_height) = Gridcell; % Initialize object array
    maxes = size(grid);

    n_cpus = 1; % Number of CPUs in system
    cpu(n_cpus) = Heatspreader(); % Initialize CPU array
    cpu_temps = zeros(n_cpus, num_points + 1);
    cpu_pos = zeros(n_cpus, 4);

    % Initialize CPU heat spreader
    cpu(1).T = ambient_ke;
    cpu(1).heat_gen = ke_inc / cpu(1).area;
    cpu(1).width = 20;
    cpu(1).height = 20;
    cpu_temps(1, 1) = cpu.T;
    cpu(1).handle = 1;
    cpu_pos(1, :) = [half_width half_height cpu(1).width cpu(1).height];


    for k = 1:length(cpu)
        % Add CPU to grid system
        for x_pos = cpu_pos(k, 1):(cpu_pos(k, 1) + cpu(k).width)
            for y_pos = cpu_pos(k, 2):(cpu_pos(k, 2) + cpu(k).height)
                grid(x_pos, y_pos).cpu_h = k;
            end
        end
    end

    wall(1) = Case_molecule(); % Empty wall list

    % Generate initial particles inside case
    min_bounds = [1 1];
    max_bounds = [round(canvas_width / 2) canvas_height];
    for i = 1:init_n_inside
        % Update particles array
        particles(i) = particles(i).rand_pos(grid, min_bounds, max_bounds);
        particles(i) = particles(i).rand_vel();
        particles(i).T = ambient_ke;
        particles(i).handle = i;
        particles(i).marker = 'square';

        % Update grid system
        grid(particles(i).rpos(1), particles(i).rpos(2)).particle_h = ...
            i;
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = ...
            false;
    end
    % Count for number of objects for which properties have been assigned
    count_assigned = init_n_inside;

    % Generate number of particles outside case
    min_bounds = [round(canvas_width / 2) 1];
    max_bounds = [canvas_width canvas_height];
    for i = (count_assigned + 1):(count_assigned + init_n_outside)
        particles(i) = particles(i).rand_pos(grid, min_bounds, max_bounds);
        particles(i) = particles(i).rand_vel();
        particles(i).handle = i;

        % Update grid system
        grid(particles(i).rpos(1), particles(i).rpos(2)).particle_h = ...
            i;
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = ...
            false;
    end
    count_assigned = count_assigned + init_n_outside; % Update counter

    % Save kinetic energy distribution of particles before
    list_KEs = zeros(1, length(particles));
    for i = 1:length(particles)
        list_KEs(i) = particles(i).T;
    end

    times = 0:dt:(end_time + dt); % Initialize array for time at each point
    % Initialize arrays for number of particles inside/outside case
    n_outside = zeros(1, num_points + 1);
    n_outside(1) = init_n_outside;
    n_inside = zeros(1, num_points + 1);
    n_inside(1) = init_n_inside;
    j = 0; % Initialize counter
    for t = times(1:num_points)
        j = j + 1; % Increment point counter
        % Print message for current time in simulation
        if (mod(round(t, 10), report_interval) == 0)
            fprintf('Time: %1.2f \n', t); % Print current time 
            fprintf('Heat gen: %1.3f KE / A \n', cpu(1).heat_gen);
        end

        n_inside(j + 1) = n_inside(j);
        n_outside(j + 1) = n_outside(j);
        for i = 1:length(particles) % Iterate through particles array
            old_x = particles(i).pos(1); % Previous x coordinate

            % Update positions of particle
            % Remove old position
            grid(particles(i).rpos(1), particles(i).rpos(2)).empty = true;

            [particles, particles(i), cpu] = ...
                particles(i).move(particles, dt, grid, cpu, maxes, wall);
            % Update new position
            grid(particles(i).rpos(1), particles(i).rpos(2)).empty = false;
            grid(particles(i).rpos(1), ...
                particles(i).rpos(2)).particle_h = ...
                particles(i).handle;

            new_x = particles(i).pos(1);
            if and(old_x > half_width, new_x <= half_width)
                n_inside(j + 1) = n_inside(j + 1) + 1;
                n_outside(j + 1) = n_outside(j + 1) - 1;
            elseif and(old_x < half_width, new_x >= half_width)
                n_inside(j + 1) = n_inside(j + 1) - 1;
                n_outside(j + 1) = n_outside(j + 1) + 1;
            end
        end

        for k = 1:length(cpu)
            % Add heat to CPU
            cpu(k) = cpu(k).add_heat(dt);
            cpu_temps(k, j + 1) = cpu(k).T;

        end
    end
    plot(times, cpu_temps(1, :), 'DisplayName', ...
        sprintf('Heat gen: %1.3f KE / A', cpu(1).heat_gen));
end
legend;

