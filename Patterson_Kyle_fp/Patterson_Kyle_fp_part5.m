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

output_vid = false; % Alternative is to render in real time
render_plots = false; % Turning this off should speed up simulation

if and(render_plots, output_vid)
    % Output AVI
    filename = strcat(part_name, '.avi');
    % Set up the movie.
    vid = VideoWriter(filename);
    vid.FrameRate = 60; % Video frame rate
    open(vid); % Open file
end

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
    
    % Draw edge
    divider = plot(ones(1, 2) * (canvas_width / 2), ...
         [0 canvas_height], 'k--');
end

n_cpus = 1; % Number of CPUs in system
cpu(n_cpus) = Heatspreader(); % Initialize CPU array
cpu_temps = zeros(n_cpus, num_points + 1);
cpu_pos = zeros(n_cpus, 4);

% Initialize CPU heat spreader
cpu(1).T = ambient_ke;

% Constant increase of overall KE on CPU per second
ke_inc = 10;
cpu(1).heat_gen = ke_inc / cpu(1).area;
cpu(1).width = 20;
cpu(1).height = 20;
cpu_temps(1, 1) = cpu.T;
cpu(1).handle = 1;
cpu_pos(1, :) = [half_width half_height cpu(1).width cpu(1).height];

% Initialize GPU heat spreader
% cpu(2).T = ambient_ke;
% ke_inc = 150; % Constant increase of overall KE on CPU per second
% cpu(2).heat_gen = ke_inc / cpu(2).area;
% cpu(2).width = 40;
% cpu(2).height = 10;
% cpu_temps(2, 1) = cpu.T;
% cpu(2).handle = 1;
% cpu_pos(2, :) = [half_width (half_height / 2) cpu(2).width cpu(2).height];

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
                (half_height/2 + cpu(k).height/2), ...
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
        
        [particles, particles(i), cpu] = ...
            particles(i).move(particles, dt, grid, cpu, maxes, wall);
        % Update new position
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = false;
        grid(particles(i).rpos(1), particles(i).rpos(2)).particle_h = ...
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
title('CPU Temperature Over Time')
xlabel('Time (seconds)');
ylabel('Temperature');

figure(); % Make a new figure for number of particles on each side
title('Gradient');
xlabel('Time (seconds)');
xlim([times(1) times(num_points + 1)]);
ylim([0 total_particles]);
ylabel('Number of particles');
hold on;
plot(times, n_inside);
plot(times, n_outside, '--');
legend('Number on left', 'Number on right');

