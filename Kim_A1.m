% ASSIGNMENT 1: Simplified, Two-Joint Arm Simulation 
% Convention: In configuration vectors, elts 1 and 2 respectively refer to  
% shoulder and elbow 

clear variables;
clc;
clf;

% Set parameters 
dur = 10; % duration of run in s
Dt = 0.01; % time step in s
psi_1 = 2.00; % state dyamics coefficient 
psi_2 = 0.63; % "
psi_3 = 0.52; % "
psi_1_est = 1*psi_1; % agent's estimate of psi_1 
psi_2_est = 1*psi_2; % " psi_2 
psi_3_est = 1*psi_3; % " psi_3 
wait = ceil(1/Dt); % interval btw. target jumps in time steps
q_min = [-2; 0]; % min. shoulder & elbow joint-angle vector in rads
q_max = [1.5; 2.5]; % max. "

% Initialise
q_star = [0; 0]; % target joint-angle vector  
q = q_star; % initial joint-angle vector 
q_vel = [0; 0]; % " velocity 
q_acc = [0; 0]; % " acceleration 
q_est = q; % model estimate of q is intially correct
q_vel_est = q_vel; % model estimate of q_vel is initially correct
hpr_1 = -13; % Hurwitz polynomial root: (x-hpr)^2 has roots hpr 
hpr_2 = -20; % " 
alpha = [2*hpr_1 -hpr_1^2; 2*hpr_2 -hpr_2^2]; % coeffs. of desired dynamics 
M = [psi_1+2*psi_2*cos(q(2)) psi_3+psi_2*cos(q(2)); ...
    psi_3+psi_2*cos(q(2)) psi_3]; % state dyamics matrix  
M_est = [psi_1_est + 2*psi_2_est*cos(q_est(2)) psi_3_est + psi_2_est*cos(q_est(2)); ... 
    psi_3_est + psi_2_est*cos(q_est(2)) psi_3_est]; % " with agent estimates
GAMMA = psi_2*sin(q(2))*[-q_vel(2) -(q_vel(1)+q_vel(2)); q_vel(1) 0]; % state dyn. matrix 
GAMMA_est = psi_2_est*sin(q_est(2))*[-q_vel_est(2) -(q_vel_est(1) + q_vel_est(2)); ... 
    q_vel_est(1) 0]; % " with agent estimates 
a = M*q_acc + GAMMA*q_vel; % action holds initial configuration 
step = 0; % # time steps since simulation began 

% Set up graphics
DATA = zeros(9, 1 + floor(dur/Dt)); % allocate memory
i_gp = 1; % index of graph pts
DATA(:, i_gp) = [0; q_star; q; q_est; a]; 

% Run time loop 
for t = Dt:Dt:dur
    % Set target 
    step = step + 1; 
    if mod(step, wait) == 1
        q_star = [q_min(1) + (q_max(1)-q_min(1))*rand(); ... 
            q_min(2) + (q_max(2)-q_min(2))*rand()]; % target joint-angles 
        q_star_dists = {-(2 + q_star(1)), 1.5-q_star(1), -q_star(2), ... 
            2.5-q_star(2)}; % distances of q_star cmpts to their bounds 
        % randomly choose one of first two and one of second two dists: 
        epsilon = [q_star_dists{randi([1, 2])}; q_star_dists{randi([3, 4])}];
        q_vel_star = [epsilon(1)*rand(); epsilon(2)*rand()]; % target velocities 
        q_est = q; % correct estimates at start of each target jump

    end 
    q_star = q_star + q_vel_star*Dt; % Euler integration 

    % Compute command using internal feedback  
    e = q_est - q_star; % joint-angle config. error estimated from model
    % implementing Hurwitz desired dynamics with sliding mode control 
    q_acc_star = alpha(:, 1).*(q_vel_est - q_vel_star) + alpha(:, 2).*e;

    %q_vel_star_target = alpha(:, 2).*e; 

    %q_acc_star = alpha(:, 1).*sign(q_vel_est - alpha(:, 2).*e);
    
    %q_acc_star = alpha(:, 1).*q_vel_est + alpha(:, 2).*e; 
    a = M_est*q_acc_star + GAMMA_est*q_vel_est;
    q_acc_est = inv(M_est)*(a - GAMMA_est*q_vel_est); % internal feedback
    q_est = q_est + Dt*q_vel_est; % Euler integration 
    q_vel_est = q_vel_est + q_acc_est*Dt; % "

    % Update state dynamics matrices with agent estimates 
    M_est = [psi_1_est + 2*psi_2_est*cos(q_est(2)) psi_3_est + psi_2_est*cos(q_est(2)); ... 
    psi_3_est + psi_2_est*cos(q_est(2)) psi_3_est]; 
    GAMMA_est = psi_2_est*sin(q_est(2))*[-q_vel_est(2) -(q_vel_est(1) + q_vel_est(2)); ... 
    q_vel_est(1) 0];  

    % Compute next state 
    q_acc = inv(M)*(a-GAMMA*q_vel); % state dynamics 
    q = q + Dt*q_vel; % Euler integration 
    q_vel = q_vel + Dt*q_acc; % " 
    q_bounded = max(q_min, min(q_max, q)); % bound joints w/in motion ranges 
    q_vel = q_vel - q_vel.*((q == q_max).*(q_vel>0)) + (q==q_min).*(q_vel < 0);
    q = q_bounded;

    % Update state dynamics matrices 
    M = [psi_1+2*psi_2*cos(q(2)) psi_3+psi_2*cos(q(2)); ...
    psi_3+psi_2*cos(q(2)) psi_3]; 
    GAMMA = psi_2*sin(q(2))*[-q_vel(2) -(q_vel(1)+q_vel(2)); q_vel(1) 0];

    % Record data for plotting 
    i_gp = i_gp + 1;
    DATA(:, i_gp) = [t; q_star; q; q_est; a]; 

    % Check action bounding 
    for i = 1:2
        if abs(a(i)) > 2000
            disp(['Component ' num2str(i) ' of action out of bounds!'])
        end 
   end 
end 

DATA = DATA(:, 1:i_gp); % discard any unused coloumns 

% Plot
figure(1); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(gcf, 'Name', 'Simplified Two-Joint Arm Simulation', 'NumberTitle', 'off'); 
subplot(2, 1, 1);  
plot(DATA(1, :), DATA(2, :), 'r:', 'LineWidth', 1);
hold on; 
grid on; 
plot(DATA(1, :), DATA(4, :), 'r'); 
plot([0, dur], [q_max(1), q_max(1)], 'color', [.5 .5 .5], 'LineStyle', ...
    '--', 'LineWidth', 0.25);
plot([0, dur], [q_min(1), q_min(1)], 'color', [.5 .5 .5], 'LineStyle', ...
    '--', 'LineWidth', 0.25);
plot(DATA(1, :), DATA(3, :), 'b:', 'LineWidth', 1); 
plot(DATA(1, :), DATA(5,:), 'b'); 
plot([0, dur], [q_max(2), q_max(2)], 'color', [.5 .5 .5], 'LineStyle', ...
    '--', 'LineWidth', 0.25);
plot([0, dur], [q_min(2), q_min(2)], 'color', [.5 .5 .5], 'LineStyle', '--', ...
    'LineWidth', 0.25);
%plot(DATA(1, :), DATA(6, :), 'r--'); 
%plot(DATA(1, :), DATA(7, :), 'b--'); 
ylim(1.05*[min(q_min), max(q_max)]); 
ylabel('q'); 
%legend('$q^{*}_{1}$', '$q_1$', '$q_1$ bds', '', '$q^{*}_{2}$', '$q_2$', ...
%    '$q_2$ bds', '', '$q_{est., 1}$', '$q_{est., 2}$');
set(gca, 'TickLength', [0, 0]); 

subplot(2, 1, 2);
plot(DATA(1, :), DATA(8, :), 'r');
hold on;
grid on; 
plot(DATA(1, :), DATA(9, :), 'b');
plot([0, dur], [2000, 2000], 'color', [.5 .5 .5], 'LineStyle', '--', ...
    'LineWidth', 0.25);
plot([0, dur], [-2000, -2000], 'color', [.5 .5 .5], 'LineStyle', '--', ...
    'LineWidth', 0.25);
xlabel('t');
ylabel('action');
ylim(1.05*[-2000, 2000]); 
%legend('$a_{1}$', '$a_{2}$'); 
set(gca, 'TickLength', [0, 0]); 

