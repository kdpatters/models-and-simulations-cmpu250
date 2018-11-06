%                        CMPU250 - Professor Eric Aaron
%                             HW2 - Kyle Patterson
%                                  April 2018

%   ####################################################################
% ###                                                                  ###
% #                               2. Brainpower                          #
% ###                                                                  ###
%   ####################################################################

% ------------------------------------------------------------------------
% Part (d)
my_title = 'Ex2d: HH, with modified Na-K switch point';

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
% Initialize arrays for concentrations
Kconce_i = zeros(1, num_points); % Concentrat. of K inside membrane (mM/L)
Kconce_o = zeros(1, num_points); % Concentrat. of K outside membrane (mM/L)
Naconce_i = zeros(1, num_points); % Concentration Na inside membrane (mM/L)
Naconce_o = zeros(1, num_points); % Concentrat. Na outside membrane (mM/L)

% Initial values
n(1) = 0.317; % Potassium activation gating variable inital value;
              % probability of K gate being open
m(1) = 0.05; % Sodium activation gating variable initial value;
             % probability of Na gate being open
h(1) = 0.6; % Sodium inactivation gating variable initial value;
            % probability of Na gate being inactivated
V(1) = -65; % Action potential initial value (mV)

% Initial concentrations
Kconce_i(1) = 150; % Init concentration of K inside membrane (mM/L)
Kconce_o(1) = 5.5; % Init concentration of K outside membrane (mM/L)
Naconce_i(1) = 15; % Init concentration Na inside membrane (mM/L)
Naconce_o(1) = 150; % Init concentration Na outside membrane (mM/L)

% Constants
% Volume calculation
eps_0 = 8.854187817e-12; % Permittivity_of_vacuum (F/m)
N_A = 6.022140857e23; % Avogadro's constant
% 10^-12 Coulombs to elementary charge conversion factor
conv_fac = 1.602176620898e7 ;
thickness = 6.3; % Average thickness of lipid bilayer (nm) 

% Calculate volumes inside/outside
conce_pos_ions_outside = Kconce_o(1) + Naconce_o(1); % mM/L
conce_pos_ions_inside = Kconce_i(1) + Naconce_i(1); % mM/L
% Additional concentration of positive ions inside
diff_conce = conce_pos_ions_inside - conce_pos_ions_outside; % mM/L
% Calculate moles of additional positive ions within cell membrane
potential = -V(1) * 10^4;
mol_pos_ions = potential * 4 * pi * eps_0 * thickness * conv_fac / N_A;
mmol_pos_ions = mol_pos_ions * 10^3; % Convert to milli moles.
vol_inside = mmol_pos_ions / diff_conce; % Volume inside cell (L)

d_squid_axon = 500e-6; % Diammeter of squid axon (m)
vol_axon_portion = d_squid_axon^3 * 1000; % Volume (L)

hp_threshold = V(1); % Hyperpolarization threshold
C_M = 0.1; % Capacitance (uF/cm^2)
V_K = -77; % Displacement from equilibrium potential for K^+ (mV)
V_L = -54.4; % Displacement from equilibrium potential for leakage (mV)
V_Na = 50; % Displacement from equilibrium potential for Na^+ (mV)
V_Na_open = -55; % Voltage at which Na gate opens (mV)
V_Na_close = 50; % Voltage at which Na gate closes (mV)
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
I_K = @(V, m, h, n, t, K_toggle) K_toggle ...
    * (g_K * n ^ 4 * (V - V_K)); 
% Leakage current (nA)
I_L = @(V, m, h, n, t, toggle) toggle * g_L * (V - V_L); 
% Sodium channel current (nA)
I_Na = @(V, m, h, n, t, Na_toggle) Na_toggle ...
    * (g_Na * m ^ 3 * h * (V - V_Na));
% Na-K pump
I_L_toggle = true;
I_P = @(toggle) toggle * I_L(V(1), m(1), h(1), n(1), start_time, ...
    I_L_toggle);

% Rate of change of `n` (ms^-1)
dndt = @(V, m, h, n, t) a_n(V) * (1 - n) - b_n(V) * n; 
% Rate of change of `m` (ms^-1)
dmdt = @(V, m, h, n, t) a_m(V) * (1 - m) - b_m(V) * m; 
% Rate of change of `h` (ms^-1)
dhdt = @(V, m, h, n, t) a_h(V) * (1 - h) - b_h(V) * h;

% Rate of change of action potential (mV/ms)
dVdt = @(V, m, h, n, t, Na_toggle, K_toggle, I_P_toggle, I_L_toggle) ...
    (I_ext(t) - I_Na(V, m, h, n, t, Na_toggle) ...
    - I_K(V, m, h, n, t, K_toggle) ...
    - I_L(V, m, h, n, t, I_L_toggle) + I_P(I_P_toggle)) / C_M;

i = 1;
prior_to_peak = true; % Open Na gate first then K gate when Na closes
Na_toggle = false; % Start with both channels closed
K_toggle = false;
I_P_toggle = true;
for time = times(2:num_points)
    i = i + 1;
    
    % Turn Na gate on
    if ((V(i - 1) >= V_Na_open) && prior_to_peak)
        Na_toggle = true; % Turn Na channel on
        K_toggle = not (Na_toggle); % Turn K channel off if not already
    end
    
    % Simulate each point using RK4 approximation
    [m(i), h(i), n(i), V(i)] = RK4(m(i - 1), dmdt, h(i - 1), ...
        dhdt, n(i - 1), dndt, V(i - 1), ...
        dVdt, time, dt, Na_toggle, K_toggle, I_P_toggle, I_L_toggle);
    
    % Concentrations fall below 0
    if (Kconce_o(i-1) <= 0)
        I_L_toggle = false;
        I_P_toggle = false;
        %K_toggle = false;
    elseif (Naconce_i(i-1) <= 0)
        I_L_toggle = true;
        I_P_toggle = false;
    else
        I_L_toggle = true;
        I_P_toggle = true;
    end
    %if (Naconce_i(i-1) <= 0)
    %    Na_toggle = false;
    %end
    
        % Change in number of K ions inside cell (C e-12)
    leakage = -I_L(V(i), m(i), h(i), n(i), time, I_L_toggle);
    K_pump = -2 * I_P(I_P_toggle); % Potassium added from Na-K pump
    %delta_K_coulombs = (leakage + K_pump + I_K(V(i), m(i), h(i), n(i), time, ...
    %    K_toggle)) * dt;
    %delta_K_moles = delta_K_coulombs * conv_fac / N_A;
    %delta_K = delta_K_moles / vol_inside;
    delta_K = (leakage + K_pump - I_K(V(i), m(i), h(i), n(i), time, ...
        K_toggle)) * dt;
    
    % Change in number of Na ions inside cell (nA * ms)
    Na_pump = 3 * I_P(I_P_toggle); % Sodium removed from Na-K pump
    %delta_Na_coulombs = (Na_pump + I_Na(V(i), m(i), h(i), n(i), time, ...
    %    Na_toggle)) * dt;
    %delta_Na_moles = delta_Na_coulombs * conv_fac / N_A;
    %delta_Na = delta_Na_moles / vol_inside;
    delta_Na = (Na_pump + I_Na(V(i), m(i), h(i), n(i), time, ...
        Na_toggle)) * dt;
    
    Kconce_o(i) = Kconce_o(i-1) - delta_K;
    Kconce_i(i) = Kconce_i(i-1) + delta_K;
    Naconce_o(i) = Naconce_o(i-1) - delta_Na;
    Naconce_i(i) = Naconce_i(i-1) + delta_Na;
    
    % Turn Na gate off and K gate on
    if ((V(i - 1) >= V_Na_close) && prior_to_peak)
        
        Na_toggle = false; % Turn Na channel off
        K_toggle = not (Na_toggle); % Turn K channel on
        prior_to_peak = false; % Voltage should decrease now
    % Have we gone below equilibrium potential?
    elseif ((V(i) >= V(i - 1)) && (V(i) < hp_threshold))
        Na_toggle = false; % Turn both channels off again
        K_toggle = false;
    end
end
                  
% Plot
my_figure = figure();

plot(times, V, '--b', 'DisplayName', 'Membrane potential (mV)');

title(my_title)
ylabel('Membrane potential (mV)')
xlabel('Time (ms)')
legend('show','Location','northeast')

gating_variables = figure();
hold on;
plot(times, m, 'r', 'DisplayName', 'Sodium activation gating variable');
plot(times, h, '-b', 'DisplayName', 'Sodium inactivation gating variable');
plot(times, n, '-.g', 'DisplayName', ...
    'Potassium activation gating variable');
title(my_title)
ylabel('Probability')
xlabel('Time (ms)')
legend('show','Location','northeast')

concentrations = figure();
hold on;
plot(times, Naconce_i, '-.b', 'DisplayName', ...
    'Intracellular Na concentration');
plot(times, Naconce_o, 'b', 'DisplayName', ...
    'Extracellular Na concentration');
plot(times, Kconce_i, '-.r', 'DisplayName', ...
    'Intracellular K concentration');
plot(times, Kconce_o, 'r', 'DisplayName', ...
    'Extracellular K concentration');
title(my_title)
ylabel('mM/L')
xlabel('Time (ms)')
legend('show','Location','northeast')
                  
% ### RK4 Approximation ###
function [m_out, h_out, n_out, V_out] = ...
    RK4(m, dmdt, h, dhdt, n, dndt, V, dVdt, t, dt, Na_toggle, ...
    K_toggle, I_P_toggle, I_L_toggle)
        t1 = t + 0.5 * dt; % Halfway through interval
        t2 = t + dt; % End of interval
        
        M1 = dmdt(V, m, h, n, t) * dt;
        H1 = dhdt(V, m, h, n, t) * dt;
        N1 = dndt(V, m, h, n, t) * dt;
        v1 = dVdt(V, m, h, n, t, Na_toggle, K_toggle, ...
            I_P_toggle, I_L_toggle) * dt;
        
        M2 = dmdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        H2 = dhdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        N2 = dndt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1) * dt;
        v2 = dVdt(V + v1 / 2, m + M1 / 2, h + H1 / 2, n + N1 / 2, t1, ... 
            Na_toggle, K_toggle, I_P_toggle, I_L_toggle) * dt;
        
        M3 = dmdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        H3 = dhdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        N3 = dndt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1) * dt;
        v3 = dVdt(V + v2 / 2, m + M2 / 2, h + H2 / 2, n + N2 / 2, t1, ...
            Na_toggle, K_toggle, I_P_toggle, I_L_toggle) * dt;
            
        M4 = dmdt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        H4 = dhdt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        N4 = dndt(V + v3, m + M3, h + H3, n + N3, t2) * dt;
        v4 = dVdt(V + v3, m + M3, h + H3, n + N3, t2, ...
            Na_toggle, K_toggle, I_P_toggle, I_L_toggle) * dt;
        
        m_out = m + (M1 + 2 * M2 + 2 * M3 + M4) / 6;
        h_out = h + (H1 + 2 * H2 + 2 * H3 + H4) / 6;
        n_out = n + (N1 + 2 * N2 + 2 * N3 + N4) / 6;
        V_out = V + (v1 + 2 * v2 + 2 * v3 + v4) / 6;
end

