%                        CMPU250 - Professor Eric Aaron
%                             HW1 - Kyle Patterson

% ------------------------------------------------------------------------
%                            1. Radioactive Decay!
% ------------------------------------------------------------------------
fprintf("\n1. Radioactive Decay!\n");
RADIUM226_DECAY = 4.27869 * 10^(-4);

% Part (a)
fprintf("\nPart (a)\n")
time_domain = 0:0.25:10^5; % in years
initial_amt = 1; % fraction of atoms remaining, starting with 100%
amts_left = rad_decay(initial_amt, RADIUM226_DECAY, time_domain);

plot(time_domain, amts_left);

% Part (b)
fprintf("\nPart (b)\n")
duration = 5*10^2; % in years
initial_amt = 1; % fraction of atoms remaining, starting with 100%
years_500 = rad_decay(initial_amt, RADIUM226_DECAY, duration);
fprintf("\nThe fraction remaining of the original quantity after" + ...  
    " %d years\nis approximately %0.5f.\n", duration, years_500);

duration = 5*10^3; % in years
years_5000 = rad_decay(initial_amt, RADIUM226_DECAY, duration);
fprintf("\nThe fraction remaining of the original quantity after" + ...  
    " %d years\nis approximately %0.5f.\n", duration, years_5000);

% Part (c)
fprintf("\nPart (c)\n")
% Our model for radioactive decay states that:
% fraction remaining == exp(-rate * time)
%
% Therefore, manipulationg this formula, we get:
% 0.60 == exp(-rate * age)
% age == log(0.60) / (-rate)
frac_remaining = 0.60;
age_60_percent_left = log(frac_remaining) / (-RADIUM226_DECAY);
fprintf("\nWhen there is the fraction %0.2f remaining of the original" + ...
    "\nquantity, approximately %0.f years have passed.\n", frac_remaining, ...
    age_60_percent_left);

% -rad_decay-
% Models radioactive decay given initial quantity, decay rate, and time
% length during which decay occurs.  Returns amount left after given
% length of time.
% :initial_amount: float representing initial amount
% :decay_rate: float representing rate at which given atoms decay.  Should
% be a number between 0 and 1 exclusive.
function amt_left = rad_decay(initial_amt, decay_rate, time)
    amt_left = initial_amt * exp(-decay_rate * time);
end