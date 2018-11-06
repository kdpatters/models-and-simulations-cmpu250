%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (6)
clear;
part_name = 'part6';

%% Parameters

% Output parameters
% Size of simulation environment (points/pixels)
canvas_width = 100;
half_width = round(canvas_width / 2);
quarter_width = (half_width / 2);
canvas_height = 100;
half_height = round(canvas_height / 2);
quarter_height = round(half_height / 2);

end_time = 20.0; % Length of simulation (seconds)
dt = 0.01; % Timestep (seconds)
num_points = ceil(end_time / dt) + 1;
% How fast time should change in the simulation relative to real time
playback_speed = 0.1;
pause_time = dt / playback_speed;

report_interval = 0.1; % Spacing between text output reports (seconds)
img_interval = dt; % Spacing between image outputs for GIF (seconds)

output_vid = false; % Alternative is to render in real time
render_plots = true; % Turning this off should speed up simulation

if and(render_plots, output_vid)
    % Output AVI
    filename = strcat(part_name, '.avi');
    % Set up the movie.
    vid = VideoWriter(filename);
    vid.FrameRate = 60; % Video frame rate
    open(vid); % Open file
end

% Ambient temperature
temp = Air_molecule;
ambient_ke = temp.ambient_ke;
max_ke = temp.max_ke;

% Initial numbers of particles
% Number of particles inside case
init_n_inside = 30;

% Number of particles outside of case
init_n_outside = 30;

total_particles = init_n_inside + init_n_outside;
particles(total_particles) = Air_molecule; % Initialize empty object array

% Create table to store positions of all particles
grid(canvas_width, canvas_height) = Gridcell; % Initialize object array
maxes = size(grid);

% Draw plot
if render_plots
    canvas = figure(); % Initialize figure
    hold on;
    axis manual;
    axis([0 canvas_width 0 canvas_height]); % Set dimensions
end

%% CPUs
n_cpus = 2; % Number of CPUs in system
cpu(n_cpus) = Heatspreader(); % Initialize CPU array
cpu_temps = zeros(n_cpus, num_points + 1);
cpu_pos = zeros(n_cpus, 4);

% Initialize CPU heat spreader
cpu(1).T = ambient_ke;
ke_inc = 10; % Constant increase of overall KE on CPU per second
cpu(1).heat_gen = ke_inc / cpu(1).area;
cpu(1).width = 20;
cpu(1).height = 20;
cpu_pos(1, :) = [half_width half_height cpu(1).width cpu(1).height];

% Initialize GPU heat spreader
cpu(2).T = ambient_ke;
ke_inc = 150; % Constant increase of overall KE on CPU per second
cpu(2).heat_gen = ke_inc / cpu(2).area;
cpu(2).width = 40;
cpu(2).height = 10;
cpu_pos(2, :) = [half_width quarter_height cpu(2).width cpu(2).height];

for i = 1:n_cpus
    cpu(i).handle = i;
    cpu_temps(i, 1) = cpu.T;
end

if render_plots
    cpu_plot(n_cpus) = rectangle('Position', cpu_pos(n_cpus, :),...
        'FaceColor', [.3 .3 .3], ...
        'EdgeColor', [.3 .3 .3], ...
        'LineWidth', 3);
end

for k = 1:length(cpu)
    if render_plots
        % Draw CPU heat spreader
        cpu_plot(k) = rectangle('Position', cpu_pos(k, :),...
            'FaceColor', [.3 .3 .3], ...
            'EdgeColor', [.3 .3 .3], ...
            'LineWidth', 3);
        if (k == 1)
            text(half_width + cpu(k).width/6, ...
                (half_height + cpu(k).height/2), ...
                'CPU', 'Color', 'white' , 'FontSize', 24);
        elseif (k == 2)
            text(half_width + cpu(k).width/6, ...
                (quarter_height + cpu(k).height/2), ...
                'GPU', 'Color', 'white' , 'FontSize', 16);
        end
    end

    % Add CPU to grid system
    for x_pos = cpu_pos(k, 1):(cpu_pos(k, 1) + cpu(k).width)
        for y_pos = cpu_pos(k, 2):(cpu_pos(k, 2) + cpu(k).height)
            grid(x_pos, y_pos).cpu_h = k;
        end
    end
end

%% Walls
n_walls = 8; % Number of walls in system
wall(n_walls) = Case_molecule(); % Initialize wall array
wall_temps = zeros(n_walls, num_points + 1);
wall_pos = zeros(n_walls, 4);

% Hardcoded walls (another way to do this would be to have a file storing
% the coordinates)

% Set positions and sizes of walls
% Coordinates of bottom left corner of case
case_pos = [(half_width - 26) (quarter_height - 10)];
case_height = 0;
case_width = 0;
wall(1).rpos = case_pos;
wall(1).width = 1;
wall(1).height = 20;
case_height = case_height + wall(1).height;

v_gap = 25;
wall(2).rpos = wall(1).rpos + [0 (wall(1).height + v_gap)];
wall(2).width = 1;
wall(2).height = 20;
case_height = case_height + v_gap + wall(2).height;

wall(3).rpos = wall(2).rpos + [0 wall(2).height];
wall(3).width = 30;
wall(3).height = 1;
case_width = case_width + wall(3).width;

h_gap = 20;
wall(4).rpos = wall(3).rpos + [(wall(3).width + h_gap) 0];
wall(4).width = 20;
wall(4).height = 1;
case_width = case_width + h_gap + wall(4).width;

wall(5).width = 1;
wall(5).height = 20;
wall(5).rpos = wall(4).rpos + [(wall(4).width - 1) -wall(5).height];

v_gap = 25;
wall(6).width = 1;
wall(6).height = 20;
wall(6).rpos = wall(5).rpos + [0 -(wall(6).height + v_gap)];

wall(7).width = 30;
wall(7).height = 1;
wall(7).rpos = wall(6).rpos + [-(wall(7).width - 1) 0];

h_gap = 20;
wall(8).width = 20;
wall(8).height = 1;
wall(8).rpos = wall(7).rpos + [-(wall(8).width + h_gap) 0]; 

% Coordinates of case sides
top_side_y = case_pos(2) + case_height;
bottom_side_y = case_pos(2);
left_side_x = case_pos(1);
right_side_x = case_pos(1) + case_width;

case_area = case_width * case_height;
outside_area = canvas_width * canvas_height - case_area;

% Whether a given (x, y) is in the PC case
inside = @(x, y) and(and(y < top_side_y, y > bottom_side_y), ...
            and(x > left_side_x, x < right_side_x));

% Whether a given (x, y) is outside the PC case        
outside = @(x, y) not(inside(x, y));
        
for i = 1:n_walls
    wall(i).therm_cond = 0; % Make the walls non-thermally conductive
    wall(i).handle = i;
    wall_temps(i, 1) = wall.T;
    wall_pos(i, :) = [wall(i).rpos wall(i).width wall(i).height];
end

if render_plots
    wall_plot(n_walls) = rectangle('Position', wall_pos(n_walls, :),...
        'FaceColor', [.3 .3 .3], ...
        'EdgeColor', [.3 .3 .3], ...
        'LineWidth', 3);
end

for k = 1:length(wall)
    if render_plots
        % Draw wall heat spreader
        wall_plot(k) = rectangle('Position', wall_pos(k, :), ...
            'FaceColor', [.3 .3 .3], ...
            'EdgeColor', [.3 .3 .3], ...
            'LineWidth', 3);
    end

    % Add wall to grid system
    for x_pos = wall_pos(k, 1):(wall_pos(k, 1) + wall(k).width)
        for y_pos = wall_pos(k, 2):(wall_pos(k, 2) + wall(k).height)
            
            grid(x_pos, y_pos).wall_h = k;
            grid(x_pos, y_pos).empty = false;
        end
    end
end

%% Generate initial particles inside case
min_bounds = case_pos;
max_bounds = case_pos + [case_width case_height];
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
for i = (count_assigned + 1):(count_assigned + init_n_outside)
    box = randi(4); % Choose one of four boxes other than in case
    if and(box == 1, case_pos(1) > 1) % Area to left of case
        min_bounds = [1 1];
        max_bounds = [(case_pos(1) - 1) canvas_height];
    elseif (box == 2) % Area above case
        min_bounds = case_pos + [0 (case_height + 1)];
        max_bounds = [(case_pos(1) + case_width) canvas_height];
    elseif (box == 3) % Area to right of case
        min_bounds = [(case_pos(1) + case_width + 1) 1];
        max_bounds = [canvas_width canvas_height];
    elseif (box == 4) % Area below case
        min_bounds = [case_pos(1) 1];
        max_bounds = [(case_pos(1) + case_width) case_pos(2)];
    else
        min_bounds = [1 1];
        max_bounds = [canvas_width canvas_height];
    end
    
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

if render_plots
    % Array of handles to plotted particles
    plot_p(total_particles) = plot(0, 0); % Initialize array
    for i = 1:total_particles % Set marker for each particle
        plot_p(i) = plot(0, 0, 'o'); % Replace placeholders with line objects
        plot_p(i).MarkerEdgeColor = [0 0 0];
    end
end

% Save kinetic energy distribution of particles before
list_KEs = zeros(1, length(particles));
for i = 1:length(particles)
    list_KEs(i) = particles(i).T;
end

%% Run simulation
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
    end
    
    n_inside(j + 1) = n_inside(j);
    n_outside(j + 1) = n_outside(j);
    
    for i = 1:length(particles) % Iterate through particles array
        old_x = particles(i).pos(1); % Previous x coordinate
        old_y = particles(i).pos(2); % Previous y coordinate
            
        if render_plots
            % Plot new positions of particles
            plot_p(i).XData = particles(i).pos(1);
            plot_p(i).YData = particles(i).pos(2);
        
            % Color new temperatures of particles
            % Basic color gradation
            temp_intensity = particles(i).T / max_ke;
            if (temp_intensity <= 1)
                p_color = [1 (1 - temp_intensity) 0];
            else
                p_color = [1 0 0]; % Max out at red
            end

            plot_p(i).MarkerFaceColor = p_color;
            plot_p(i).Marker = particles(i).marker;
        end
        
        % Update positions of particle
        % Remove old position
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = true;
        
        [particles, particles(i), cpu, wall] = ...
            particles(i).move(particles, dt, grid, cpu, maxes, wall);
        % Update new position
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = false;
        grid(particles(i).rpos(1), particles(i).rpos(2)).particle_h = ...
            particles(i).handle;
        
        new_x = particles(i).pos(1);
        new_y = particles(i).pos(2);
        
        % Particle leaves the case
        if and(inside(old_x, old_y), outside(new_x, new_y))
            n_inside(j + 1) = n_inside(j + 1) - 1;
            n_outside(j + 1) = n_outside(j + 1) + 1;
        % Particle enters the case
        elseif and(inside(new_x, new_y), outside(old_x, old_y))
            n_inside(j + 1) = n_inside(j + 1) + 1;
            n_outside(j + 1) = n_outside(j + 1) - 1;
        end
    end
    
    for k = 1:length(cpu)
        % Add heat to CPU
        cpu(k) = cpu(k).add_heat(dt);
        cpu_temps(k, j + 1) = cpu(k).T;

        % Redraw CPU
        if render_plots
            cpu_temp_intensity = cpu(k).T / max_ke;
            if (cpu_temp_intensity <= 1)
                cpu_color = [1 (1 - cpu_temp_intensity) 0];
            else
                cpu_color = [1 0 0]; % Max out at red
            end
            cpu_plot(k).EdgeColor = cpu_color;
            cpu_plot(k).FaceColor = cpu_color;
        end
    end
    
%    for k = 1:length(wall)
%         wall_temps(k, j + 1) = wall(k).T;
%        
%         % Redraw walls
%         if render_plots
%             wall_temp_intensity = wall(k).T / max_ke;
%             if (wall_temp_intensity <= 1)
%                 wall_color = [1 (1 - wall_temp_intensity) 0];
%             else
%                 wall_color = [1 0 0]; % Max out at red
%             end
%             wall_plot(k).EdgeColor = wall_color;
%             wall_plot(k).FaceColor = wall_color;
%         end
%     end
    
    if and(render_plots, output_vid)
        drawnow; % Call draw function
        
        if (mod(t, img_interval) == 0)
            frame = getframe(gcf);
            writeVideo(vid, frame);
        end
    elseif render_plots
        pause(pause_time); % Call draw function then wait
    end
end
if and(render_plots, output_vid)
    close(vid); % Stop writing to AVI file
end

%% Draw final plots
% Make histogram for distribution of temperatures
ke_dist_init = figure();
histogram(list_KEs, 'BinWidth', 5);
title('Distribution of KEs before simulation')
xlabel('Kinetic energy');
ylabel('Number');

ke_dist_new = figure();
% Make a new figure for kinetic energy distribution of particles
list_KEs_new = zeros(1, length(particles));
for i = 1:length(particles)
    list_KEs_new(i) = particles(i).T;
end
histogram(list_KEs_new, 'BinWidth', 5);
title('Distribution of KEs after simulation')
xlabel('Kinetic energy');
ylabel('Number');

% Make a new figure for kinetic energy distribution of particles
cpu_heat = figure();
hold on;
for k = 1:length(cpu)
    plot(times, cpu_temps(k, :));
end
legend('CPU', 'GPU');
title('CPU/GPU Temperature Over Time')
xlabel('Time (seconds)');
ylabel('Temperature');

figure(); % Make a new figure for number of particles on each side
title('Pressure Inside/Outside Case');
xlabel('Time (seconds)');
xlim([times(1) times(num_points + 1)]);
ylabel('Number of particles per unit area');
hold on;
plot(times, n_inside / case_area);
plot(times, n_outside / outside_area, '--');
legend('Pressure inside', 'Pressure outside');

