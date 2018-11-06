%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                            1. Have a Ball!                           #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (c)

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

% Initial conditions
positions(1) = 11; % Initial position above ground (m)
velocities(1) = 15; % Initial velocity (m/s)

% Forces
weight = mass * g; % Force of gravity on ball (N)
total_force = @(p, v, t) weight;

% Acceleration
acc = @(p, v, t) g + 0.01 * (v + p) + 0.3 * t ^ 2;
% Initial acceleration (m/s^2)
accelerations(1) = acc(positions(1), velocities(1), 0);
% Velocity
vel = @(p, v, a, t) (a - g - 0.01 * p - 0.3 * t ^ 2) / 0.01;

i = 1;
for time = times(2:num_points)
    i = i + 1;
    [positions(i), velocities(i)] = RK4(positions(i - 1), ... 
        vel, velocities(i - 1), acc, time, dt);
    accelerations(i) = acc(positions(i), velocities(i), time);
end

% Plot
my_figure = figure();
hold on;

yyaxis left;
ylabel('Position (m)')
plot(times, positions, '--b', 'DisplayName', 'Position (m)');

yyaxis right;
ylabel('Velocity (m/s) and Acceleration (m/s^2)')
plot(times, velocities, 'r', 'DisplayName', 'Velocity (m/s)');

plot(times, accelerations, '-.g', 'DisplayName', 'Acceleration (m/s^2)');

title('Ex1c: Vertical throw using RK4')
xlabel('Time (seconds)')
legend('show','Location','southeast')

% ### RK4 Approximation ###
function [p_out, v_out] = ...
    RK4(pt, v, vt, a, t, dt)
        t1 = t + 0.5 * dt; % Halfway through interval
        t2 = t + dt; % End of interval
        
        A1 = a(pt, vt, t);
        A1dt = A1 * dt;
        V1dt = v(pt, vt, A1, t) * dt;
        
        A2 = a(pt + V1dt / 2, vt + A1dt / 2, t1);
        A2dt = A2 * dt;
        V2dt = v(pt + V1dt / 2, vt + A1dt / 2, A2, t1) * dt;
        
        A3 = a(pt + V2dt / 2, vt + A2dt / 2, t2);
        A3dt = A3 * dt;
        V3dt = v(pt + V2dt / 2, vt + A2dt / 2, A3, t2) * dt;
            
        A4 = a(pt + V3dt, vt + A3dt, t2);
        A4dt = A4 * dt;
        V4dt = v(pt + V3dt, vt + A3dt, A4, t2) * dt;
        
        p_out = pt + (V1dt + 2 * V2dt + 2 * V3dt + V4dt) / 6;
        v_out = vt + (A1dt + 2 * A2dt + 2 * A3dt + A4dt) / 6;
end








