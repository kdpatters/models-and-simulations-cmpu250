%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                            1. Change It Up!                        #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (a)

% ### Parameters of simulation ###
dt = 0.01; % Units in months
months_sim = 12; % Number of months to run simulation over
num_points = months_sim / dt; % Number of points to simulate over

% Declare simulation arrays
Y = zeros(1, num_points);
P = zeros(1, num_points);
H = zeros(1, num_points);

% Population of Tuna
Y(1) = 1000; % Initial population in thousands of tuna
Y_grow = 2; % Tuna growth rate coefficient
Y_prey_P = 2e-2; % Tuna death to predation coefficient

% Population of Sharks
P(1) = 150; % Initial population
P_death = 1.06; % Shark death rate coefficient
P_feed_Y = 1e-2; % Shark growth due to feeding on tuna
P_fish = 5e-2; % Shark death due to humans

% Population of Humans
H(1) = 10;
H_death = 1; % Human death due to population size
H_feed_P = 1e-2; % Human population growth from feeding on sharks

% Change in tuna population
dYdt = @(Yt, Pt, Ht, t) Y_grow * Yt - Y_prey_P * (Pt * Yt);

% Change in shark population
dPdt = @(Yt, Pt, Ht, t) - P_death * Pt + P_feed_Y * (Pt * Yt) ...
       - P_fish * (Ht * Pt);

% Change in human population
dHdt = @(Yt, Pt, Ht, t) - H_death * Ht + H_feed_P * (Ht * Pt);   

points = 1:num_points;
times = points * dt;

% Simulation loop
for j = points(1:num_points - 1)
    [Y(j + 1), P(j + 1), H(j + 1)] = ...
        RK4(dYdt, Y, dPdt, P, dHdt, H, times, dt, j);
end

% Plot
populations = figure();
hold on;

% Plot tuna
plot(times, Y, 'b', 'DisplayName', 'Tuna');

% Plot sharks
plot(times, P, '--r', 'DisplayName', 'Sharks');

% Plot sharks
plot(times, H, '-.g', 'DisplayName', 'Humans');

title(strcat('Human death =', int2str(H_death)))
xlabel('Time (times)')
ylabel('Population')
legend('show')

H_death = 0.5;

% Change in human population
dHdt = @(Yt, Pt, Ht, t) - H_death * Ht + H_feed_P * (Ht * Pt); 

% Simulation loop
for j = points(1:num_points - 1)
    [Y(j + 1), P(j + 1), H(j + 1)] = ...
        RK4(dYdt, Y, dPdt, P, dHdt, H, times, dt, j);
end

%Plot
populations = figure();
hold on;

% Plot tuna
plot(times, Y, 'b', 'DisplayName', 'Tuna');

% Plot sharks
plot(times, P, '--r', 'DisplayName', 'Sharks');

% Plot sharks
plot(times, H, '-.g', 'DisplayName', 'Humans');

title(strcat('Human death =', num2str(H_death)))
xlabel('Time (times)')
ylabel('Population')
legend('show')

% ### RK4 Approximation ###
function [Y_out, P_out, H_out] = ...
    RK4(dYdt, Y, dPdt, P, dHdt, H, times, dt, j)
        t0 = times(j); % Start of interval
        t1 = (j + 0.5) * dt; % Halfway through interval
        t2 = times(j + 1); % End of interval
        
        A1 = dYdt(Y(j), P(j), H(j), t0) * dt;
        B1 = dPdt(Y(j), P(j), H(j), t0) * dt;
        C1 = dHdt(Y(j), P(j), H(j), t0) * dt;
        
        A2 = dYdt(Y(j) + A1 / 2, P(j) + B1 / 2, H(j) + C1 / 2, t1) * dt;
        B2 = dPdt(Y(j) + A1 / 2, P(j) + B1 / 2, H(j) + C1 / 2, t1) * dt;
        C2 = dHdt(Y(j) + A1 / 2, P(j) + B1 / 2, H(j) + C1 / 2, t1) * dt;
        
        A3 = dYdt(Y(j) + A2 / 2, P(j) + B2 / 2, H(j) + C2 / 2, t1) * dt;
        B3 = dPdt(Y(j) + A2 / 2, P(j) + B2 / 2, H(j) + C2 / 2, t1) * dt;
        C3 = dHdt(Y(j) + A2 / 2, P(j) + B2 / 2, H(j) + C2 / 2, t1) * dt;
        
        A4 = dYdt(Y(j) + A3 / 2, P(j) + B3 / 2, H(j) + C3 / 2, t2) * dt;
        B4 = dPdt(Y(j) + A3 / 2, P(j) + B3 / 2, H(j) + C3 / 2, t2) * dt;
        C4 = dHdt(Y(j) + A3 / 2, P(j) + B3 / 2, H(j) + C3 / 2, t2) * dt;
        
        Y_out = Y(j) + (A1 + 2 * A2 + 2 * A3 + A4) / 6;
        P_out = P(j) + (B1 + 2 * B2 + 2 * B3 + B4) / 6;
        H_out = H(j) + (C1 + 2 * C2 + 2 * C3 + C4) / 6;
end

