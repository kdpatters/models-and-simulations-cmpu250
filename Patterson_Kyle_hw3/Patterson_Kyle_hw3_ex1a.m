%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                            1. Have a Ball!                           #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (a)

% Constants
% Assign constants; away from earth is positive direction
mass = 0.5; % Mass (kg)
g = -9.81; % Acceleration due to gravity (m/s^2)
end_time = 4; % Time in seconds
start_time = 0; % Time in seconds
dt = 0.25; % Timestep in seconds
% Number points to simulate
num_points = ceil((end_time - start_time) / dt) + 1;

% Assign empty arrays of zeros to store velocities and positions at
% each point
velocities = zeros(1, num_points);
positions = zeros(1, num_points);
accelerations = zeros(1, num_points);
times = start_time:dt:end_time;

% Initial velocity and position
velocities(1) = 15; % Initial velocity (m/s)
positions(1) = 11; % Initial position above ground (m)

% Forces
weight = mass * g; % Force of gravity on ball (N)
total_force = @(p, v, a, t) weight;

acc = @(p, v, a, t) total_force(p, v, a, t) / mass; % Constant acceleration
accelerations(1) = acc(positions(1), ...
    velocities(1), accelerations(1), 0); % Initial acceleration (m/s^2)
vel = @(p, v, a, t) v + a * dt; % Linear velocity
pos = @(p, v, a, t) p + v * dt;
true_pos = @(p, v, a, t) positions(1) + velocities(1)*t + ...
    accelerations(1) / 2 * t.^2;

i = 1;
for time = times(2:num_points)
    i = i + 1;
    accelerations(i) = acc(positions(i - 1), velocities(i - 1), ...
        accelerations(i - 1), time);
    velocities(i) = vel(positions(i - 1), velocities(i - 1), ...
        accelerations(i - 1), time);
    positions(i) = pos(positions(i - 1), velocities(i - 1), ...
        accelerations(i - 1), time);
end

% Plot
my_figure = figure();
hold on;

plot(times, positions, '--b', 'DisplayName', 'Position (m)');
plot(times, true_pos(positions, velocities, accelerations, ...
    times), '-.g', 'DisplayName', 'True position (m)');

plot(times, velocities, 'r', 'DisplayName', 'Velocity (m/s)');

title("Ex1a: Vertical throw using Euler's method")
xlabel('Time (seconds)')
legend('show')







