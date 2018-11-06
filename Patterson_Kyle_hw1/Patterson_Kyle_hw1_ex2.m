%                        CMPU250 - Professor Eric Aaron
%                             HW1 - Kyle Patterson

% ------------------------------------------------------------------------
%                            2. A Big Piece of Pi!
% ------------------------------------------------------------------------

% Part (a)
fprintf("\nPart (a)\n")
num_points = 1000;
estimate_1000 = approx_pi(num_points);
fprintf("\nWith %d simulations, we get an estimate of %0.6f for pi.\n", ...
     num_points, estimate_1000);
 
% Create a logarithmically spaced array of length `num_points` 
% from `min_sims` to `max_sims`
min_sims = 0; % 10^0, so this will be a minimum of 1
max_sims = 6;
num_points = 10;
sizes_of_sims = logspace(min_sims, max_sims, num_points);

num_sims = 10^2; % We'll do 100 simulations at each size
deviations = zeros([1,num_points]); % Initialize array to hold stddev at each size
means = zeros([1,num_points]); % Initialize array to hold sample means
curr_estimates = zeros([1,num_points]);
for size_ind = 1:num_points
    for k = 1:num_sims
        curr_estimates(k) = approx_pi(sizes_of_sims(size_ind));
    end
    deviations(size_ind) = std(curr_estimates);
    means(size_ind) = mean(curr_estimates);
end
figure();
loglog(sizes_of_sims, deviations);
ylim([0 max(deviations)])
title('Effect of number of points on stddev for approx of Pi')
xlabel('Number of points')
ylabel('Standard deviation')

figure();
pi_y(1:sizes_of_sims) = pi;
semilogx(sizes_of_sims, means, sizes_of_sims, pi_y);
title('Effect of number of points on mean for approx Pi')
xlabel('Number of points')
ylabel('Sample mean')

% Part (b)
fprintf("\nPart (b)\n")

most_dense_plot = 4; % 10^4
least_dense_plot = 1; % 10^1
max_plots = 5;
for trials = round(logspace(max([min_sims least_dense_plot]), ... 
    min([max_sims most_dense_plot]), ... 
    min([num_points max_plots])))
    x = rand([1, trials]); % Generate x coordinates
    y = rand([1, trials]); % Generate y coordinates

    figure();
    hold on;
    for i = 1:trials
        if in_circle(x(i), y(i))
            plot(x(i), y(i), '.r')
        else
            plot(x(i), y(i), '.b')
        end
    end
    title(['Approximation with ', num2str(trials), ' simulations'])
end

% Part (c)
fprintf("\nPart (c)\n")

tests = 10^3;
chosen_n = 10^4;
values_pi = zeros([1,tests]);

for i = 1:tests
    values_pi(i) = approx_pi(chosen_n);
end

maximum = max(values_pi);
minimum = min(values_pi);
average = mean(values_pi);

fprintf("\nWith %d tests and a chosen number of simulations of %d," + ...
    "\nwe get a min of %0.5f, a max of %0.5f, and an average of" + ...
    "\n%0.5f", tests, chosen_n, minimum, maximum, average);

% Part (d)
fprintf("\nPart (d)\n")

% -approx-
% Given number, runs that number of Monte Carlo simulations to approximate
% value of pi.  Returns float representing approximate value of pi.
% :n: integer number of simulations
function approx = approx_pi(n)
    n = round(n); % Ensure `n` is integer
    x = rand([2, n]); % Generate x coordinates
    y = rand([2, n]); % Generate y coordinates
    count_inside = 0;
    
    for i = 1:n
        count_inside = count_inside + in_circle(x(i), y(i));
    end
    num_quadrants = 4; % Find area of whole circle
    approx = num_quadrants * (count_inside / n);
end

% -inside-
% Given x and y coordinates, returns bool for whether given point is either
% both positive y quadrant of unit circle.
% :x: float representing x coordinate
% :y: float representing y coordinate
function inside = in_circle(x, y)
    inside = (y <= unit_circle_y(x));
end

% -unit_circle_y-
% Given x coordinate, returns corresponding positive y coordinate of unit
% circle.
% :x: float representing x coordinate
function y = unit_circle_y(x)
    y = sqrt(1 - x^2);
end

