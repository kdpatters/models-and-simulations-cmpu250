%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                           1. Up The Food Chain!                      #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (b), (c)

% ### Parameters of simulation ###
dt = 0.01; % Units in months
months_sim = 12; % Number of months to run simulation over
num_points = months_sim / dt; % Number of points to simulate over

% Declare simulation arrays
Y = zeros(1, num_points);
P = zeros(1, num_points);

% Population of Tuna
Y(1) = 100; % Initial population in thousands of tuna
Y_grow = 2; % Tuna growth rate coefficient
Y_prey_P = 2e-2; % Tuna death to predation coefficient
Y_fish = 1; % Tuna death due to human coefficient

% Population of Sharks
P(1) = 15; % Initial population
P_death = 1.06; % Shark death rate coefficient
P_feed_Y = 1e-2; % Shark growth due to feeding on tuna
P_fish = Y_fish; % Shark death due to humans

% Change in tuna population
dYdt = @(Yt, Pt, Ht, t) Y_grow * Yt - Y_prey_P * (Pt * Yt) ...
       - Y_fish * (Ht * Yt);

% Change in shark population
dPdt = @(Yt, Pt, Ht, t) - P_death * Pt + P_feed_Y * (Pt * Yt) ...
       - P_fish * (Ht * Pt);
   
% Change in human population
dHdt = @(Yt, Pt, Ht, t) 0;   

% ### Effect from humans ###
% Fishing rate params
min_fr = 0.01;
max_fr = 2;

H = @(t) 0;

num_fr = 3;

% Growth rates to test for tuna
min_gr_t = 0.01;
max_gr_t = 1;
num_gr_t = 3;

% Growth rates to test for sharks
min_gr_s = 0.01;
max_gr_s = 0.1;
num_gr_s = 3;

% Death rates to test for sharks and tuna
min_dr = 0.01;
max_dr = 1;
num_dr = 3;

points = 1:num_points;
times = points * dt;

%H = 0; % Ignore human component for first few simulations
   
% Try different growth rates for tuna
myInterval = (max_gr_t - min_gr_t)/(num_gr_t - 1);
for gr = min_gr_t:myInterval:max_gr_t
    
    dYdt2 = @(Yt, Pt, Ht, t) gr * Yt - Y_prey_P * (Pt * Yt);

    % Simulation loop
    for j = points(1:num_points - 1)
        [Y(j + 1), P(j + 1)] = ...
            RK4(dYdt2, Y, dPdt, P, H, times, dt, j);
    end

    % Plot
    populations = figure();
    hold on;

    % Plot tuna
    plot(times, Y, 'b', 'DisplayName', 'Tuna');

    % Plot sharks
    plot(times, P, '--r', 'DisplayName', 'Sharks');

    title(strcat('Tuna growth rate=', ' ', num2str(gr)))
    xlabel('Time (times)')
    ylabel('Population')
    legend('show')
end

% Try different growth rates for shark
myInterval = (max_gr_s - min_gr_s)/(num_gr_s - 1);
for gr = min_gr_s:myInterval:max_gr_s
    
    dPdt2 = @(Yt, Pt, Ht, t) - P_death * Pt + gr * (Pt * Yt);

    % Simulation loop
    for j = points(1:num_points - 1)
        [Y(j + 1), P(j + 1)] = ...
            RK4(dYdt, Y, dPdt2, P, H, times, dt, j);
    end

    % Plot
    populations = figure();
    hold on;

    % Plot tuna
    plot(times, Y, 'b', 'DisplayName', 'Tuna');

    % Plot sharks
    plot(times, P, '--r', 'DisplayName', 'Sharks');

    title(strcat('Shark growth rate=', ' ', num2str(gr)))
    xlabel('Time (times)')
    ylabel('Population')
    legend('show')
end

% Try different death rates for tuna
for dr = min_dr:(max_dr - min_dr)/(num_dr - 1):max_dr
    
    dYdt3 = @(Yt, Pt, Ht, t) Y_grow * Yt - dr * (Pt * Yt);

    % Simulation loop
    for j = points(1:num_points - 1)
        [Y(j + 1), P(j + 1)] = ...
            RK4(dYdt3, Y, dPdt, P, H, times, dt, j);
    end

    % Plot
    populations = figure();
    hold on;

    % Plot tuna
    plot(times, Y, 'b', 'DisplayName', 'Tuna');

    % Plot sharks
    plot(times, P, '--r', 'DisplayName', 'Sharks');

    title(strcat('Tuna death rate=', ' ', num2str(dr)))
    xlabel('Time (months)')
    ylabel('Population')
    legend('show')
end

% Try different death rates for shark
for dr = min_dr:(max_dr - min_dr)/(num_dr - 1):max_dr
    
    dPdt3 = @(Yt, Pt, Ht, t) - dr * Pt + P_feed_Y * (Pt * Yt);

    % Simulation loop
    for j = points(1:num_points - 1)
        [Y(j + 1), P(j + 1)] = ...
            RK4(dYdt, Y, dPdt3, P, H, times, dt, j);
    end

    % Plot
    populations = figure();
    hold on;

    % Plot tuna
    plot(times, Y, 'b', 'DisplayName', 'Tuna');

    % Plot sharks
    plot(times, P, '--r', 'DisplayName', 'Sharks');

    title(strcat('Shark death rate=', ' ', num2str(dr)))
    xlabel('Time (months)')
    ylabel('Population')
    legend('show')
end

% Try different fishing rates
for fr = min_fr:(max_fr - min_fr)/(num_fr - 1):max_fr % Test fishing rates
    
    % Set new fishing rate
    H = @(t) fr;
    
    % Simulation loop
    for j = points(1:num_points - 1)
        [Y(j + 1), P(j + 1)] = ...
            RK4(dYdt, Y, dPdt, P, H, times, dt, j);
    end
    
    % Plot
    populations = figure();
    hold on;
    
    % Plot tuna
    plot(times, Y, 'b', 'DisplayName', 'Tuna');

    % Plot sharks
    plot(times, P, '--r', 'DisplayName', 'Sharks');
    
    title(strcat('Fishing rate=', ' ', num2str(fr)))
    xlabel('Time (months)')
    ylabel('Population')
    legend('show')
end

% ### RK4 Approximation ###
function [Y_out, P_out] = RK4(dYdt, Y, dPdt, P, H, times, dt, j)
        t0 = times(j); % Start of interval
        t1 = (j + 0.5) * dt; % Halfway through interval
        t2 = times(j + 1); % End of interval
        
        A1 = dYdt(Y(j), P(j), H(t0), t0) * dt;
        B1 = dPdt(Y(j), P(j), H(t0), t0) * dt;
        
        A2 = dYdt(Y(j) + A1 / 2, P(j) + B1 / 2, H(t1), t1) * dt;
        B2 = dPdt(Y(j) + A1 / 2, P(j) + B1 / 2, H(t1), t1) * dt;
        
        A3 = dYdt(Y(j) + A2 / 2, P(j) + B2 / 2, H(t1), t1) * dt;
        B3 = dPdt(Y(j) + A2 / 2, P(j) + B2 / 2, H(t1), t1) * dt;
       
        A4 = dYdt(Y(j) + A3 / 2, P(j) + B3 / 2, H(t2), t2) * dt;
        B4 = dPdt(Y(j) + A3 / 2, P(j) + B3 / 2, H(t2), t2) * dt;
        
        Y_out = Y(j) + (A1 + 2 * A2 + 2 * A3 + A4) / 6;
        P_out = P(j) + (B1 + 2 * B2 + 2 * B3 + B4) / 6;
end

