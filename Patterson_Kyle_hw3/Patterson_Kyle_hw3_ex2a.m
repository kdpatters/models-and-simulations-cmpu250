%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                               2. Brainpower                          #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (a)

% Experimental parameters
dt = 0.001; % Timestep (ms)
start_time = 0; % Time at start of simulation (ms)
end_time = 3; % Time at end of simulation (ms)
% Number points to simulate
num_points = ceil((end_time - start_time) / dt) + 1;

% Initialize arrays for gating variables
n = zeros(1, num_points); % Potassium activation gating variable
m = zeros(1, num_points); % Sodium activation gating variable
h = zeros(1, num_points); % Sodium inactivation gating variable
% Initialize array from voltages
V = zeros(1, num_points); % Action potential (mV)
% Initialize array for time at each point evaluated
times = start_time:dt:end_time;

% Initial values
n(1) = 0.317; % Potassium activation gating variable inital value;
              % probability of K gate being open
m(1) = 0.05; % Sodium activation gating variable initial value;
             % probability of Na gate being open
h(1) = 0.6; % Sodium inactivation gating variable initial value;
            % probability of Na gate being inactivated
V(1) = -65; % Action potential initial value (mV)

% Constants
hp_threshold = V(1); % Hyperpolarization threshold
C_M = 0.1; % Capacitance (uF/cm^2)
V_K = -77; % Displacement from equilibrium potential for K^+ (mV)
V_L = -54.4; % Displacement from equilibrium potential for leakage (mV)
V_Na = 50; % Displacement from equilibrium potential for Na^+ (mV)
V_Na_open = -55; % Voltage at which Na gate opens (mV)
V_Na_close = 49.3; % Voltage at which Na gate closes (mV)
V_K_open = V_Na_close; % Voltage at which K gate opens (mV)
g_K = 35; % Maximum K conductance (mS/cm^2)
g_Na = 120; % Maximum Na conductance (mS/cm^2)
g_L = 0.3; % Maximum leakage conductance (mS/cm^2)

% Variables
% Opening rate constant (ms^-1)
a_n = @(V) 0.01 * (V + 55) / (1 - exp(-(V + 55) / 10));
% Opening rate constant (ms^-1)
a_m = @(V) 0.1 * (V + 40) / (1 - exp(-(V + 40) / 10));
a_h = @(V) 0.07 * exp(-(V + 65) / 20); % Opening rate constant (ms^-1)
b_n = @(V) 0.125 * exp(-(V + 65) / 20); % Closing rate constant (ms^-1)
b_m = @(V) 4 * exp(-(V + 65) / 80); % Closing rate constant (ms^-1)
b_h = @(V) 1 / (exp(-(V + 35) / 10) + 1); % Closing rate constant (ms^-1)

stimulus_start = 0.5; % Start of stimulus current (ms)
stimulus_duration = 0.5; % Length of stimulus current (ms)
stimulus_end = stimulus_start + stimulus_duration;
% Stimulus current (mV); piecewise function
I_ext = @(t) (stimulus_start <= t) * (t <= stimulus_end) * 15;
% Potassium channel current (nA)
I_K = @(V, m, h, n, t) g_K * n ^ 4 * (V - V_K); 
I_L = @(V, m, h, n, t) g_L * (V - V_L); % Leakage current (nA)
% Sodium channel current (nA)
I_Na = @(V, m, h, n, t) (g_Na * m ^ 3 * h * (V - V_Na));

% Rate of change of `n` (ms^-1)
dndt = @(V, m, h, n, t) a_n(V) * (1 - n) - b_n(V) * n; 
% Rate of change of `m` (ms^-1)
dmdt = @(V, m, h, n, t) a_m(V) * (1 - m) - b_m(V) * m; 
% Rate of change of `h` (ms^-1)
dhdt = @(V, m, h, n, t) a_h(V) * (1 - h) - b_h(V) * h;

% Rate of change of action potential (mV/ms)
dVdt = @(V, m, h, n, t) (I_ext(t) - I_L(V, m, h, n, t)) / C_M;

i = 1;
Na_on_K_off = true; % Open Na gate first then K gate when Na closes
for time = times(2:num_points)
    i = i + 1;
    
    % Turn Na gate on
    if ((V(i - 1) >= V_Na_open) && Na_on_K_off)

        % Rate of change of action potential (mV/ms)
        dVdt = @(V, m, h, n, t) (I_ext(t) - I_Na(V, m, h, n, t) ...
            - I_L(V, m, h, n, t)) / C_M;
    end
    
    % Simulate each point using RK4 approximation
    [m(i), h(i), n(i), V(i)] = RK4(m(i - 1), dmdt, h(i - 1), ...
        dhdt, n(i - 1), dndt, V(i - 1), ...
        dVdt, time, dt);
    
    % Turn Na gate off and K gate on
    if ((V(i - 1) >= V_Na_close) && Na_on_K_off)
        
         % Rate of change of action potential (mV/ms)
        dVdt = @(V, m, h, n, t) (I_ext(t) - I_K(V, m, h, n, t) ...
                      - I_L(V, m, h, n, t)) / C_M;
        Na_on_K_off = false;   
     % Have we gone below equilibrium potential?
    elseif ((V(i) >= V(i - 1)) && (V(i) < hp_threshold))
        dVdt = @(V, m, h, n, t) (I_ext(t) - I_L(V, m, h, n, t)) / C_M;
    end
end
                  
% Plot
my_figure = figure();

plot(times, V, '--b', 'DisplayName', 'Membrane potential (mV)');

title('Ex2a: HH, no ATPase pump')
ylabel('Membrane potential (mV)')
xlabel('Time (ms)')
legend('show','Location','northeast')

gating_variables = figure();
hold on;
plot(times, m, 'r', 'DisplayName', 'Sodium activation gating variable');
plot(times, h, '-b', 'DisplayName', 'Sodium inactivation gating variable');
plot(times, n, '-.g', 'DisplayName', ...
    'Potassium activation gating variable');
title('Ex2a: HH, no ATPase pump')
ylabel('Probability')
xlabel('Time (ms)')
legend('show','Location','northeast')
                  
% ### RK4 Approximation ###
function [m_out, h_out, n_out, V_out] = ...
    RK4(m, dmdt, h, dhdt, n, dndt, V, dVdt, t, dt)
        t1 = t + 0.5 * dt; % Halfway through interval
        t2 = t + dt; % End of interval
        
        M1 = dmdt(V, m, h, n, t) * dt;
        H1 = dhdt(V, m, h, n, t) * dt;
        N1 = dndt(V, m, h, n, t) * dt;
        v1 = dVdt(V, m, h, n, t) * dt;
        
        M2 = dmdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        H2 = dhdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        N2 = dndt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        v2 = dVdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        
        M3 = dmdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        H3 = dhdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        N3 = dndt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        v3 = dVdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
            
        M4 = dmdt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        H4 = dhdt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        N4 = dndt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        v4 = dVdt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        
        m_out = m + (M1 + 2 * M2 + 2 * M3 + M4) / 6;
        h_out = h + (H1 + 2 * H2 + 2 * H3 + H4) / 6;
        n_out = n + (N1 + 2 * N2 + 2 * N3 + N4) / 6;
        V_out = V + (v1 + 2 * v2 + 2 * v3 + v4) / 6;
end

