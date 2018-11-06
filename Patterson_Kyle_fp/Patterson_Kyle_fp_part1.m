%                        CMPU250 - Professor Eric Aaron
%                        Final Project - Kyle Patterson
%                                  May 2018

%   ####################################################################
% ###                                                                  ###
% #                         PC Case Heat Transfer                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (1)
clear;
part_name = 'part1';

%%% Parameters

% Output parameters
% Size of simulation environment (points/pixels)
canvas_width = 100;
canvas_height = 100;

% Spacing between text output reports
report_interval = 0.1;

output_gif = false; % Alternative is to render in real time

% Output GIF filename
filename = strcat(part_name, '.gif');

end_time = 0.70; % Length of simulation (seconds)
dt = 0.01; % Timestep (seconds)
num_points = ceil(end_time / dt) + 1;
% How fast time should change in the simulation relative to real time
playback_speed = 0.1;
pause_time = dt / playback_speed;

% Ambient temperature
temp = Particle_part123();
ambient_ke = temp.ambient_ke;

% Initial numbers of particles
% Number of particles inside case
init_n_inside = 1;

% Number of particles outside of case
init_n_outside = 1;

total_particles = init_n_inside + init_n_outside;
particles(total_particles) = Particle_part123; % Initialize empty object array

% Create table to store positions of all particles
grid(canvas_width, canvas_height) = Gridcell; % Initialize object array

% Generate two particles, heading toward each other
% Generate initial particles inside case
for i = 1:init_n_inside
    % Update particles array
    particles(i).pos = [round(canvas_width / 4) canvas_height / 2];
    particles(i).rpos = round(particles(i).pos);
    particles(i).vel = [1 0] * particles(i).ke;
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
for i = (count_assigned + 1):(count_assigned + init_n_outside)
    particles(i).pos = [round(3 * canvas_width / 4) canvas_height / 2];
    particles(i).rpos = round(particles(i).pos);
    particles(i).vel = [-1 0] * particles(i).ke;
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
    
for t = 0:dt:end_time
    % Print message for current time in simulation
    if (mod(t, report_interval) == 0)
        fprintf('Time: %1.2f \n', t); % Print current time 
    end
    
    for i = 1:length(particles) % Iterate through particles array
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
    end
    
    if output_gif
        drawnow; % Call draw function
        
        % Capture current plot
        frame = getframe(canvas); 
        img = frame2im(frame); 
        [img_ind,color_map] = rgb2ind(img ,256); 

        % Write image to file
        if (t == 0) % First image
            imwrite(img_ind, color_map, ...
                filename, 'gif', 'Loopcount', inf, 'DelayTime', 0); 
        else % Additional images
            imwrite(img_ind, color_map, ...
                filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0); 
        end

    else
        pause(pause_time); % Call draw function then wait
    end
end
close; % Close figure
