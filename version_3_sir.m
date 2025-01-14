% MATLAB Code with Fixes

clc;
clear all;

%%% Exercise 1 Part A - SIR Model Simulation

% Parameters
% Infection rate
beta = 2;   
% Recovery rate
gamma = 0.1;  
% Mortality rate
mu = 0.01;      

% Initial conditions [S; I; R]
initial_conditions = [100; 1; 0]; 

% Time span (200 days)
T = [0 200];

% Solve the system using ode45
[t, X] = ode45(@(t, X) sir_func(t, X, beta, gamma, mu), T, initial_conditions);

% Extract S, I, R from the solution
S = X(:, 1);
I = X(:, 2);
R = X(:, 3);

% Plot the results
figure;
plot(t, S, 'b', t, I, 'r', t, R, 'g', 'LineWidth', 1.5);
xlabel('Time (Days)');
ylabel('Number of Individuals');
legend('Susceptible (S)', 'Infected (I)', 'Recovered (R)');
title('Deterministic SIR Model Simulation');
grid on;

% Function definition for SIR model
function dXdt = sir_func(~, X, beta, gamma, mu)
    S = X(1);
    I = X(2);
    R = X(3);
	% Total population
    N = S + I + R; 
	
    % Susceptible differential eq
    dSdt = -beta * (I * S) / N; 
	% Infected differential eq
    dIdt = beta * (I * S) / N - gamma * I - mu * I; 
	% Recovered differential eq
    dRdt = gamma * I; 
	% Putting results of all ODEs in a vector
    dXdt = [dSdt; dIdt; dRdt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Exercise 1 Part B - SIR Model Deaths Simulation

% Time span and initial population
N0 = sum(initial_conditions);

% Different beta values
figure;
hold on;
beta_values = [0.2, 0.5, 1, 2];
for i = 1:length(beta_values)
    beta_local = beta_values(i);
    [t, X] = ode45(@(t, X) sir_func(t, X, beta_local, gamma, mu), T, initial_conditions);
    N = sum(X, 2);
    deaths = N0 - N;
    plot(t, deaths, 'LineWidth', 1.5);
    disp(['Number of deaths after 200 days for beta = ', num2str(beta_local), ': ', num2str(deaths(end))]);
end

xlabel('Time (Days)');
ylabel('Number of Deaths');
legend('Beta = 0.2', 'Beta = 0.5', 'Beta = 1', 'Beta = 2');
title('Deterministic SIR Model Deaths Simulation');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Exercise 2 Part A - Stochastic SIR Model Simulations

% Two independent examples - (pass upper limit of T)
for k = 1:2 
    stochastic_sir_simulation(beta, gamma, mu, T(2), k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Exercise 2 Part B - Analysis over 1000 Samples

% Number of samples
s = 1000; 

% Part i: Mean epidemic duration
epidemic_duration(s, T(2), beta, gamma, mu); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Exercise 3 - Modified Stochastic SIR Model

% Part A: Stochastic SIR with immunity vs. susceptibility
for k = 1:2
    modified_stochastic_sir_simulation(beta, gamma, mu, T(2), k); 
end

% Part B: Repeating all Exercise 2 calculations with modifications
% Part i: Mean epidemic duration
modified_epidemic_duration(s, T(2), beta, gamma, mu);

% Part ii: Mean total number of deaths
modified_mean_deaths_calculation(s, T(2), beta, gamma, mu);

% Part iii: Mean time to epidemic peak
modified_mean_peak_time_calculation(s, T(2), beta, gamma, mu);

% Part iv: Probability of peak infected > 10
modified_probability_peak_infected(s, T(2), beta, gamma, mu);

% Part v: Vaccination rate analysis
modified_vaccination_rate_analysis(s, T(2), beta, gamma, mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Supporting Functions

function stochastic_sir_simulation(beta, gamma, mu, Tmax, k)
    S = 100;
    I = 1;
    R = 0;
    N = S + I + R;

    times = 0;
    S_values = S;
    I_values = I;
    R_values = R;

    while times(end) < Tmax && I > 0
        rate_infection = beta * S * I / N;
        rate_recovery = gamma * I;
        rate_death = mu * I;
        total_rate = rate_infection + rate_recovery + rate_death;

        tau = exprnd(1 / max(total_rate, 1e-8)); % Safeguard for low rates
        times(end + 1) = times(end) + tau;

        event_prob = rand * total_rate;

        if event_prob < rate_infection
            S = S - 1;
            I = I + 1;
        elseif event_prob < rate_infection + rate_death
            I = I - 1;
        else
            I = I - 1;
            R = R + 1;
        end

        S_values(end + 1) = S;
        I_values(end + 1) = I;
        R_values(end + 1) = R;
    end

    figure;
    plot(times, S_values, 'b', times, I_values, 'r', times, R_values, 'g', 'LineWidth', 1.5);
    xlabel('Time (Days)');
    ylabel('Population');
    legend('Susceptible (S)', 'Infected (I)', 'Recovered (R)');
    title(['Stochastic SIR Model (Trial ', num2str(k), ')']);
    grid on;
end

function epidemic_duration(s, Tmax, beta, gamma, mu)
    durations = zeros(1, s);
    for i = 1:s
        S = 100;
        I = 1;
        R = 0;
        N = S + I + R;

        times = 0;

        while I > 0 && times(end) < Tmax
            rate_infection = beta * S * I / N;
            rate_recovery = gamma * I;
            rate_death = mu * I;
            total_rate = rate_infection + rate_recovery + rate_death;

            % Safeguard for low rates
            tau = exprnd(1 / max(total_rate, 1e-8)); 
            times(end + 1) = times(end) + tau;

            event_prob = rand * total_rate;

            if event_prob < rate_infection
                S = S - 1;
                I = I + 1;
            elseif event_prob < rate_infection + rate_death
                I = I - 1;
            else
                I = I - 1;
                R = R + 1;
            end
        end

        durations(i) = times(end);
    end

    % Calculate and display mean duration
    mean_duration = mean(durations);
    disp(['Mean epidemic duration: ', num2str(mean_duration), ' days']);

    % Plot histogram
    figure;
    histogram(durations, 20, 'FaceColor', 'blue');
    xlabel('Epidemic Duration (Days)');
    ylabel('Frequency');
    title('Histogram of Epidemic Duration (Stochastic SIR Model)');
    grid on;
end


%%% Supporting Functions for Modified Stochastic SIR Model

function modified_stochastic_sir_simulation(beta, gamma, mu, Tmax, k)
    % Modified simulation with immunity considerations
    S = 100;
    I = 1;
    R = 0;
    N = S + I + R;

    times = 0;
    S_values = S;
    I_values = I;
    R_values = R;

    while times(end) < Tmax && I > 0
        rate_infection = beta * S * I / N;
        rate_recovery = gamma * I;
        rate_death = mu * I;
        total_rate = rate_infection + rate_recovery + rate_death;

        tau = exprnd(1 / max(total_rate, 1e-8)); 
        times(end + 1) = times(end) + tau;

        event_prob = rand * total_rate;

        if event_prob < rate_infection
            S = S - 1;
            I = I + 1;
        elseif event_prob < rate_infection + rate_death
            I = I - 1;
        else
            I = I - 1;
            R = R + 1;
        end

        S_values(end + 1) = S;
        I_values(end + 1) = I;
        R_values(end + 1) = R;
    end

    figure;
    plot(times, S_values, 'b', times, I_values, 'r', times, R_values, 'g', 'LineWidth', 1.5);
    xlabel('Time (Days)');
    ylabel('Population');
    legend('Susceptible (S)', 'Infected (I)', 'Recovered (R)');
    title(['Modified Stochastic SIR Model (Trial ', num2str(k), ')']);
    grid on;
end

function modified_epidemic_duration(s, Tmax, beta, gamma, mu)
    durations = zeros(1, s);
    for i = 1:s
        S = 100;
        I = 1;
        R = 0;
        N = S + I + R;

        times = 0;

        while I > 0 && times(end) < Tmax
            rate_infection = beta * S * I / N;
            rate_recovery = gamma * I;
            rate_death = mu * I;
            total_rate = rate_infection + rate_recovery + rate_death;

            tau = exprnd(1 / max(total_rate, 1e-8)); % Safeguard for low rates
            times(end + 1) = times(end) + tau;

            event_prob = rand * total_rate;

            if event_prob < rate_infection
                S = S - 1;
                I = I + 1;
            elseif event_prob < rate_infection + rate_death
                I = I - 1;
            else
                I = I - 1;
                R = R + 1;
            end
        end

        durations(i) = times(end);
    end

    % Calculate and display mean duration
    mean_duration = mean(durations);
    disp(['Modified Mean epidemic duration: ', num2str(mean_duration), ' days']);

    % Plot histogram
    figure;
    histogram(durations, 20, 'FaceColor', 'green');
    xlabel('Epidemic Duration (Days)');
    ylabel('Frequency');
    title('Histogram of Epidemic Duration (Modified SIR Model)');
    grid on;
end


function modified_mean_deaths_calculation(s, Tmax, beta, gamma, mu)
    deaths = zeros(1, s);
    for i = 1:s
        S = 100;
        I = 1;
        R = 0;
        N = S + I + R;

        while I > 0 && sum([S, I, R]) > 0
            rate_infection = beta * S * I / N;
            rate_recovery = gamma * I;
            rate_death = mu * I;
            total_rate = rate_infection + rate_recovery + rate_death;

            tau = exprnd(1 / max(total_rate, 1e-8)); 

            event_prob = rand * total_rate;

            if event_prob < rate_infection
                S = S - 1;
                I = I + 1;
            elseif event_prob < rate_infection + rate_death
                I = I - 1;
            else
                I = I - 1;
                R = R + 1;
            end
        end

        deaths(i) = 100 - (S + R);
    end

    mean_deaths = mean(deaths);
    disp(['Modified Mean total number of deaths: ', num2str(mean_deaths)]);
end

function modified_mean_peak_time_calculation(s, Tmax, beta, gamma, mu)
    peak_times = zeros(1, s);
    for i = 1:s
        S = 100;
        I = 1;
        R = 0;
        N = S + I + R;

        times = 0;
        infected_counts = I;

        while I > 0 && times(end) < Tmax
            rate_infection = beta * S * I / N;
            rate_recovery = gamma * I;
            rate_death = mu * I;
            total_rate = rate_infection + rate_recovery + rate_death;

            tau = exprnd(1 / max(total_rate, 1e-8)); 
            times(end + 1) = times(end) + tau;

            event_prob = rand * total_rate;

            if event_prob < rate_infection
                S = S - 1;
                I = I + 1;
            elseif event_prob < rate_infection + rate_death
                I = I - 1;
            else
                I = I - 1;
                R = R + 1;
            end

            infected_counts(end + 1) = I;
        end

        [~, peak_index] = max(infected_counts);
        peak_times(i) = times(peak_index);
    end

    mean_peak_time = mean(peak_times);
    disp(['Modified Mean time to epidemic peak: ', num2str(mean_peak_time), ' days']);
end

function modified_probability_peak_infected(s, Tmax, beta, gamma, mu)
    peaks_above_threshold = zeros(1, s);
    threshold = 10;
    for i = 1:s
        S = 100;
        I = 1;
        R = 0;
        N = S + I + R;

        infected_counts = I;

        while I > 0
            rate_infection = beta * S * I / N;
            rate_recovery = gamma * I;
            rate_death = mu * I;
            total_rate = rate_infection + rate_recovery + rate_death;

            tau = exprnd(1 / max(total_rate, 1e-8)); 

            event_prob = rand * total_rate;

            if event_prob < rate_infection
                S = S - 1;
                I = I + 1;
            elseif event_prob < rate_infection + rate_death
                I = I - 1;
            else
                I = I - 1;
                R = R + 1;
            end

            infected_counts(end + 1) = I;
        end

        if max(infected_counts) > threshold
            peaks_above_threshold(i) = 1;
        end
    end

    probability = mean(peaks_above_threshold);
    disp(['Modified Probability of peak infected > ', num2str(threshold), ': ', num2str(probability)]);
end

function modified_vaccination_rate_analysis(s, Tmax, beta, gamma, mu)
    vaccination_rates = 0:0.01:1;
    probabilities = zeros(1, length(vaccination_rates));

    for v = 1:length(vaccination_rates)
        vax_rate = vaccination_rates(v);
        peaks_above_threshold = zeros(1, s);
        threshold = 10;

        for i = 1:s
            S = 100 * (1 - vax_rate);
            I = 1;
            R = 0;
            N = S + I + R;

            infected_counts = I;

            while I > 0
                rate_infection = beta * S * I / N;
                rate_recovery = gamma * I;
                rate_death = mu * I;
                total_rate = rate_infection + rate_recovery + rate_death;

                tau = exprnd(1 / max(total_rate, 1e-8)); 

                event_prob = rand * total_rate;

                if event_prob < rate_infection
                    S = S - 1;
                    I = I + 1;
                elseif event_prob < rate_infection + rate_death
                    I = I - 1;
                else
                    I = I - 1;
                    R = R + 1;
                end

                infected_counts(end + 1) = I;
            end

            if max(infected_counts) > threshold
                peaks_above_threshold(i) = 1;
            end
        end

        probabilities(v) = mean(peaks_above_threshold);
    end

    base_probability = probabilities(1);
    target_probability = base_probability / 2;
    idx = find(probabilities <= target_probability, 1);

    if ~isempty(idx)
        required_vaccination_rate = vaccination_rates(idx) * 100;
        disp(['Modified Vaccination rate required to halve the peak probability: ', num2str(required_vaccination_rate), '%']);
    else
        disp('No vaccination rate found to halve the peak probability within the range.');
    end
end