%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (4)
clear;
part_name = 'part4';

%%% Parameters

% Output parameters
% Size of simulation environment (points/pixels)
canvas_width = 100;
half_width = canvas_width / 2;
canvas_height = 100;

output_gif = false; % Alternative is to render in real time

% Output GIF filename
filename = strcat(part_name, '.gif');

end_time = 2.0; % Length of simulation (seconds)
dt = 0.01; % Timestep (seconds)
report_interval = 0.1; % Spacing between text output reports (seconds)
img_interval = dt; % Spacing between image outputs for GIF (seconds)
num_points = ceil(end_time / dt) + 1;
% How fast time should change in the simulation relative to real time
playback_speed = 1;
pause_time = dt / playback_speed;

% Ambient temperature
temp = Particle;
ambient_ke = temp.ambient_ke;

% Initial numbers of particles
% Number of particles inside case
init_n_inside = 30;

% Number of particles outside of case
init_n_outside = 30;

total_particles = init_n_inside + init_n_outside;
% Initialize empty object array
particles(total_particles) = Particle_part123;

% Create table to store positions of all particles
grid(canvas_width, canvas_height) = Gridcell; % Initialize object array

% Generate initial particles inside case
min_bounds = [1 1];
max_bounds = [round(canvas_width / 2) canvas_height];
for i = 1:init_n_inside
    % Update particles array
    particles(i) = particles(i).rand_pos(grid, min_bounds, max_bounds);
    particles(i) = particles(i).rand_vel();
    particles(i).ke = ambient_ke / 10;
    particles(i).handle = i;
    particles(i).marker = 'x';
    
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

% Draw plot
canvas = figure(); % Initialize figure
hold on;
axis manual;
axis([0 canvas_width 0 canvas_height]); % Set dimensions

% Draw edge
divider = plot(ones(1, 2) * (canvas_width / 2), ...
     [0 canvas_height], 'k--');

% Array of handles to plotted particles
plot_p(total_particles) = plot(0, 0); % Initialize array
for i = 1:total_particles % Set marker for each particle
    plot_p(i) = plot(0, 0, 'ro'); % Replace placeholders with line objects
end

% Make a new figure for kinetic energy distribution of particles
list_KEs = zeros(1, length(particles));
for i = 1:length(particles)
    list_KEs(i) = particles(i).ke;
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
            
        % Plot new positions of particles
        plot_p(i).XData = particles(i).pos(1);
        plot_p(i).YData = particles(i).pos(2);
        
        % Color new temperatures of particles
        % Basic color gradation
        temp_intensity = particles(i).ke / ambient_ke;
        if (temp_intensity <= 1)
            p_color = [1 (1 - temp_intensity) 0];
        else
            p_color = [1 0 0]; % Max out at red
        end
        
        plot_p(i).MarkerEdgeColor = p_color;
        plot_p(i).MarkerFaceColor = p_color;
        plot_p(i).Marker = particles(i).marker;
        
        % Update positions of particle
        % Remove old position
        grid(particles(i).rpos(1), particles(i).rpos(2)).empty = true;
        
        [particles, particles(i)] = ...
            particles(i).move(particles, dt, grid);
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
    
    if and(output_gif, mod(t, img_interval) == 0)
        drawnow; % Call draw function
        
        % Capture current plot
        frame = getframe(canvas); 
        img = frame2im(frame); 
        [img_ind,color_map] = rgb2ind(img ,256); 

        % Write image to file
        if (t == 0) % First image
            imwrite(img_ind, color_map, ...
                filename, 'gif', 'Loopcount', inf, 'DelayTime', ...
            img_interval); 
        else % Additional images
            imwrite(img_ind, color_map, ...
                filename, 'gif', 'WriteMode', 'append', 'DelayTime', ...
                img_interval); 
        end

    else
        pause(pause_time); % Call draw function then wait
    end
end

ke_dist_new = figure();
% Make a new figure for kinetic energy distribution of particles
list_KEs_new = zeros(1, length(particles));
for i = 1:length(particles)
    list_KEs_new(i) = particles(i).ke;
end
histogram(list_KEs_new, 'BinWidth', 5);
title('Distribution of KEs after simulation')
xlabel('Kinetic energy');
ylabel('Number');

% Make histogram for distribution of temperatures
ke_dist_init = figure();
histogram(list_KEs, 'BinWidth', 5);
title('Distribution of KEs before simulation')
xlabel('Kinetic energy');
ylabel('Number');
set(0, 'currentfigure', canvas); % Make canvas active again

figure(); % Make a new figure for number of particles on each side
title('Passive potential reduction of a gradient');
xlabel('Time (seconds)');
xlim([times(1) times(num_points + 1)]);
ylim([0 total_particles]);
ylabel('Number of particles');
hold on;
plot(times, n_inside);
plot(times, n_outside, '--');
legend('Number on left', 'Number on right');

