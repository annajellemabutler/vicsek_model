%% Basic Vicsek Model Simulation
% Written by Anna Jellema-Butler
% On 12/9/2023 at UCSB 

% =========================================================
% Aim: To reproduce the Vicsek model in MATLAB and analyse the order
% parameter at different noise levels. 
% Steps of algorithm:
    % 1. Initialization
    % 2. Simulation loop
    % 3. Order parameter calculations 
    %   3.i. Plot OP over time
    %   3.ii. Plot steady-state OP for combinations of eta & N
% ==========================================================

%% 1. Initialization of parameters

% Set parameters for simulation
sim_time = 200; % max time to simulate over
vel = 1; % initial velocity
N = 1500; % particle number
L = 7 ; % size of box
r = 1; % interacting radius
eta = 0.8; % noise parameter (determines max value for random noise)
dt = 1; % time step
es = 1; % flag that indicates whether to exit loop when system stable 
prs_size = 100; % for convergence calc
prs_num = 0.01; % for convergence calc

% Initialize domain and particles with random positions and angles
x = L*rand(1,N);
y = L*rand(1,N);
theta = 2*pi*(rand(1,N)-0.5);

% Create arrays to store info
v = vel*ones(1,N); % magnitude of velocity for each particle
oP = zeros(1,sim_time); % order parameter at each time step

%% 2. Main simulation loop

% Setting up video file 
%videoObj = VideoWriter('Vicsek_Model_Simulation_condition_X.mp4', 'MPEG-4'); 
%videoObj.Quality = 100; 
%videoObj.FrameRate = 10;
%open(videoObj);

for time = 1:sim_time
    
    % Plot current particle positions, colored by angle
    % -------------------------------------------------
    scatter(x, y, 10,atan2(sin(theta), cos(theta)),'filled');
    colormap(spring);
    colorbar;
    clim([-pi, pi]);
    axis([0, L, 0, L]);
    set(gcf, 'Color', 'w')
    title(['Vicsek Model, Time Step: ' num2str(time)]);
    xlabel('X-axis');
    ylabel('Y-axis');

    % Write current frame to the video
    %writeVideo(videoObj, getframe(gcf));
    %pause(0.001);

    % Set up necesary arrays for next steps
    temp_x = zeros(1,N);
    temp_y = zeros(1,N);
    av_theta = zeros(1,N);
    phi = zeros(1,N);

    % Calculate average angle in the interacting circle
    % -------------------------------------------------
    % A) Get all pairwise distances
    dist = pdist([x' y'], 'euclidean');
    
    % B) Account for periodic boundaries 
    % - Necessary because for e.g. particles at x = 0.5 and x = 7 are
    % actually neighbors. 
    % - Thus, temporarily add L to any co-ords < r & subtract L from any
    % co-ords > L-r, and re-calculate distances. 
    % - Finally, concatenate D with temp distances and select the
    % smallest.
    temp_x(x<r) = L + x(x<r);
    temp_x(x>L-r) = x(x>L-r)-L;
    temp_x(r<=x & x<=L-r) = x(r<=x & x<=L-r);
    temp_y(y<r) = L + y(y<r);
    temp_y(y>L-r) = y(y>L-r)-L;
    temp_y(r<=y & y<=L-r) = y(r<=y & y<=L-r);
    temp_D = pdist([temp_x' temp_y'],'euclidean');   
    dist = min([dist; temp_D]); 

    % C) Convert to matrix representation
    M = squareform(dist); 

    % D) Find particles where distance < r & compute their average
    % direction angle 
    [row, col]=find(0<=M & M<r); 
    for i = 1:N
        list = row(col==i); % get list of neighbors
        av_theta(i) = atan2(mean(sin(theta(list))),mean(cos(theta(list))));
    end
        
    % Update positions, adjusting for periodic boundary conditions
    % -----------------------------------------------------------
    x = x + v.*cos(theta).*dt;
    y = y + v.*sin(theta).*dt;

    x(x<0) = L + x(x<0);
    x(L<x) = x(L<x) - L;
    y(y<0) = L + y(y<0);
    y(L<y) = y(L<y) - L;
    
    % Update angle of movement according to neighbors
    % -------------------------------------------------
    theta = av_theta + eta'.*(rand(1,N) - 0.5); 
    
    % Av. velocity calculations
    % -------------------------------------------------
    % A) Get the sum of all horizontal and verticle components of velocity 
    vel_x = sum(v.*cos(theta)); 
    vel_y = sum(v.*sin(theta)); 

    % B) Calculate the order parameter as sum of all particle velocities 
    % divided by number of velocities 
    oP(time) = ((vel_x.^2) + (vel_y^2)).^0.5 ./ sum(v);

    % C) Perform a convergence check
    % - Has the system reached a stable state?
    % - If the order parameter exhibits minimal variation over the specified
    % time window, the loop is exited. 
    if time>=prs_size % checking if enough time steps have passed 
        conv_Phi = oP(time-prs_size+1:time); % extract subset of order parameters 
        % from the time window specified
        if (max(conv_Phi)-min(conv_Phi)<prs_num) && (es==1), break,end 
    end
end

%close(videoObj);
 
if es==1 % Give v_a as the average of oP over the stable period ('steady
    % state order parameter')
    v_a = mean(conv_Phi);
    conv_ctr = abs(normalize(conv_Phi, 'center')); % subtracts the mean
    I = find(conv_ctr==min(conv_ctr),1); % finds minimum value index
    T = time-prs_size+I; % time at which minimum normalized OP occurred
else % give v_a as Phi at last time step
    v_a = oP(end)  
    T = NaN; % there is no specific stability time 
end


%% 3.i. Plotting the order parameter over time 

figure;
plot(1:sim_time, oP, 'LineWidth', 1.5)
title(['Order Parameter over Time (N=' num2str(N) ', L=' num2str(L) ', \eta=' num2str(eta) ')'])
xlabel('Frames');
ylabel('Average Velocity')
set(gcf, 'Color', 'w')
print('Va_plot_over_time', '-dpng', '-r600');

%% 3.ii. Plotting steady-state order parameter as a function of combinations of eta & N
% Specify the range of eta values to explore
eta_values = 0:0.5:5.0; % Adjust the range and step size as needed

% Define different combinations of N and L
param_combinations = [
    40, 3.1;
    100, 5;
    400, 10
];

% Initialize a cell array to store results for each combination
results_cell = cell(size(param_combinations, 1), 1);

% Loop over different combinations of N and L
for param_index = 1:size(param_combinations, 1)
    N = param_combinations(param_index, 1);
    L = param_combinations(param_index, 2);

    % Initialize arrays to store steady-state average velocities for each eta
    steady_state_velocities = zeros(size(eta_values));

    % Loop over different values of eta
    for eta_index = 1:length(eta_values)
        % Set parameters for simulation
        sim_time = 200;
        vel = 1;
        r = 1;
        eta = eta_values(eta_index);
        dt = 1;
        es = 1;
        prs_size = 100;
        prs_num = 0.01;

        % Initialize domain and particles with random positions and angles
        x = L * rand(1, N);
        y = L * rand(1, N);
        theta = 2 * pi * (rand(1, N) - 0.5);

        % Create arrays to store info
        v = vel * ones(1, N);
        Phi = zeros(1, sim_time);

        particle_trajectory = zeros(sim_time, N, 2);

        % Main simulation loop
        for time = 1:sim_time
           
        
            % Set up necesary arrays
            temp_x = zeros(1,N);
            temp_y = zeros(1,N);
            av_theta = zeros(1,N);
            phi = zeros(1,N);
        
            % Calculate average angle in the interacting circle
            % -------------------------------------------------
            % A) Get all pairwise distances
            dist = pdist([x' y'], 'euclidean');
            
            % B) Account for periodic boundaries 
            % - Necessary because for e.g. particles at x = 0.5 and x = 7 are
            % actually neighbors. 
            % - Thus, temporarily add L to any co-ords < r & subtract L from any
            % co-ords > L-r, and re-calculate distances. 
            % - Finally, concatenate D with temp distances and select the
            % smallest.
            temp_x(x<r) = L + x(x<r);
            temp_x(x>L-r) = x(x>L-r)-L;
            temp_x(r<=x & x<=L-r) = x(r<=x & x<=L-r);
            temp_y(y<r) = L + y(y<r);
            temp_y(y>L-r) = y(y>L-r)-L;
            temp_y(r<=y & y<=L-r) = y(r<=y & y<=L-r);
            temp_D = pdist([temp_x' temp_y'],'euclidean');   
            dist = min([dist; temp_D]); 
        
            % C) Convert to matrix representation
            M = squareform(dist); 
        
            % D) Find particles where distance < r & compute their average
            % direction angle 
            [row, col]=find(0<=M & M<r); 
            for i = 1:N
                list = row(col==i); % get list of neighbors
                av_theta(i) = atan2(mean(sin(theta(list))),mean(cos(theta(list))));
            end
                
            % Update positions, adjusting for periodic boundary conditions
            % -----------------------------------------------------------
            x = x + v.*cos(theta).*dt;
            y = y + v.*sin(theta).*dt;
        
            x(x<0) = L + x(x<0);
            x(L<x) = x(L<x) - L;
            y(y<0) = L + y(y<0);
            y(L<y) = y(L<y) - L;
            
            % Update angle of movement according to neighbors
            % -------------------------------------------------
            theta = av_theta + eta'.*(rand(1,N) - 0.5);
        
            % Update particle_trajectory array
            particle_trajectory(time, :, 1) = x; % store x positions
            particle_trajectory(time, :, 2) = y; % store y positions
        
            % Order parameter calculations
            % -------------------------------------------------
            % A) Get the sum of all horizontal and verticle components of velocity 
            vel_x = sum(v.*cos(theta)); 
            vel_y = sum(v.*sin(theta)); 
        
            % B) Calculate the order parameter as sum of all particle velocities 
            % divided by number of velocities 
        
            Phi(time) = ((vel_x.^2) + (vel_y^2)).^0.5 ./ sum(v);
        
            % C) Perform a convergence check
            % - Has the system reached a stable state?
            % - If the order parameter exhibits minimal variation over the specified
            % time window, the loop is exited. 
            if time>=prs_size % checking if enough time steps have passed 
                conv_Phi = Phi(time-prs_size+1:time); % extract subset of Phi 
                % from the time window specified
                if (max(conv_Phi)-min(conv_Phi)<prs_num) && (es==1), break,end 
            end
        end
    % Calculate steady-state average velocity as the average order parameter
        % over the stable period
        if es == 1
            v_a = mean(Phi(time - prs_size + 1:time));
        else
            v_a = Phi(end); % use the last value of Phi if stability is not achieved
        end

        % Store the steady-state average velocity for the current eta
        steady_state_velocities(eta_index) = v_a;
    end
    % Store the results for the current combination of N and L
    results_cell{param_index} = struct('N', N, 'L', L, 'eta_values', eta_values, 'velocities', steady_state_velocities);
end

% Plot the results for each combination of N and L
figure;
hold on;
for param_index = 1:size(param_combinations, 1)
    plot(results_cell{param_index}.eta_values, results_cell{param_index}.velocities, 'LineWidth', 1.5);
end
hold off;

xlabel('\eta', 'FontSize', 13);
ylabel('\vee_{\alpha}', 'FontSize', 13);
legend('N = 40, L = 3.1', 'N = 100, L = 5', 'N = 400, L = 10');
set(gcf, 'Color', 'w');
print('Steady_State_Velocity_vs_Eta_Multiple', '-dpng', '-r600');
