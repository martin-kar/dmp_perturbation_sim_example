%% Preliminaries
clc; close all; clear all
dmpParams;

% Simulate demonstrated trajectory:
dt = 1/250;
t = linspace(0,10,10/dt)';
traj = max(0,-sin(2*pi*t(1:round(2*end/5))/2.5)); traj = [traj ; traj(end)*ones(50,1)];

% Determine DMP from demonstration:
T_end = length(traj)*dt;
tau =T_end/3; % Yields 95% convergence at the deomnstration time
P = length(traj); % Number of time steps
w = traj2w(traj,dt, tau, c, D, alpha_z, beta_z, alpha_x, n_kernel); % Basis function weights. Used to determine f.
g = traj(end); % DMP goal point.

%% Simulate unperturbed DMP
%Initialize variables
x = 1; % Phase parameter
ydot = 0;
y = zeros(size(traj));
y(1) = traj(1);
y0 = y(1);
yddot_unpert = zeros(size(traj));

% Simulate over time t
for t = 2:2*P
    yddot  = dmp2acc(y0,y(t-1),ydot,g,tau, w, x, c, D, alpha_z, beta_z); % Get acceleration given DMP and states
    ydot = ydot + yddot*dt;
    y(t) = y(t-1) + ydot*dt;
    xdot = -alpha_x*x/tau;    
    x = x + xdot*dt;
    yddot_unpert(t) = yddot;
end
y_unpert = y;


%% Simulate DMP execution with stopping perturbation
%Initialize variables
x = 1; % Phase variable
y = zeros(size(traj));
y(1) = traj(1); % Coupled trajectory
y0 = y(1); % Start position
ya = y; % Actual trajectory
ya_dot = 0;
e = 0; % Low-pass filtered error state
ya_ddot_log = zeros(size(y));
y_ddot_log = zeros(size(y));

% Simulate over time t
for t = 2:2*P
    tau_adapt = tau*(1+(kc*e^2)); % Adaptive time parameter
    ya_ddot = get_ya_ddot_lowgain_ff(ya(t-1), ya_dot, y(t-1), ydot, yddot); % Get reference accceleration for actual trajectory
    [ydot, yddot]  = dmp2vel_acc_ss(y0, y(t-1), g, tau_adapt, w, x, dt, alpha_e, c, D, alpha_z, beta_z,ya(t-1),e,kc,tau); % Get coupled acceleration and velocit,y given DMP and states
    y(t) = y(t-1) + ydot*dt;
    xdot = -alpha_x*x/tau_adapt;    
    x = x + xdot*dt;
    
    % Stopping perturbation
    ya_dot = ya_dot + ya_ddot*dt;
    if (t > 500 && t < 750)
        ya_dot = 0;
    end
    
    ya(t) = ya(t-1) + ya_dot*dt;
    e_dot = alpha_e*(ya(t)-y(t)-e);
    e = e + e_dot*dt;
    ya_ddot_log(t) = ya_ddot;
    y_ddot_log(t) = yddot;
end
t = cumsum(dt*ones(size(y)));

% Plot results
figure
subplot(211)
plot(t(1:10:end),ya(1:10:end),'b-',t(1:10:end),y(1:10:end),'r--', t(1:10:end),y_unpert(1:10:end),'k-.','LineWidth',2)
legend('y_a','y_c','y_u')
axis([0 8 -.2 1])
xlabel('Time [s]')
ylabel('Position [m]')

subplot(212)
plot(t(1:10:end),ya_ddot_log(1:10:end),'b-',t(1:10:end),y_ddot_log(1:10:end),'r--', t(1:10:end),yddot_unpert(1:10:end),'k-.','LineWidth',2)
legend('ÿ_a','ÿ_c','ÿ_u')
axis([0 8 -20 20])
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')

%% Simulate DMP execution with moving perturbation
% Initialize variables
x = 1; % Phase variable
ydot = 0;
yddot = 0;
y = zeros(size(traj));
y(1) = traj(1);
y0 = y(1);
ya = y;
ya_dot = 0;
e = 0; % Low-pass filtered error state
ya_ddot_log = zeros(size(y));
y_ddot_log = zeros(size(y));

for t = 2:2*P
    tau_adapt = tau*(1+(kc*e^2)); % Adaptive time parameter
    ya_ddot = get_ya_ddot_lowgain_ff(ya(t-1), ya_dot, y(t-1), ydot, yddot); % Get reference accceleration for actual trajectory
    [ydot, yddot]  = dmp2vel_acc_ss(y0, y(t-1), g, tau_adapt, w, x, dt, alpha_e, c, D, alpha_z, beta_z,ya(t-1),e,kc,tau); % Get coupled acceleration and velocit,y given DMP and states
    y(t) = y(t-1) + ydot*dt;
    xdot = -alpha_x*x/tau_adapt;    
    x = x + xdot*dt;
    
    % Moving perturbation
    ya_dot = ya_dot + ya_ddot*dt;
    if (t > 500 && t < 750)
        ya_dot = .15;
    end
    
    ya(t) = ya(t-1) + ya_dot*dt;
    e_dot = alpha_e*(ya(t)-y(t)-e);
    e = e + e_dot*dt;
    ya_ddot_log(t) = ya_ddot;
    y_ddot_log(t) = yddot;
end
t = cumsum(dt*ones(size(y)));

% Plot results
figure
subplot(211)
plot(t(1:10:end),ya(1:10:end),'b-',t(1:10:end),y(1:10:end),'r--', t(1:10:end),y_unpert(1:10:end),'k-.','LineWidth',2)
legend('y_a','y_c','y_u')
axis([0 8 -.2 1.5])
xlabel('Time [s]')
ylabel('Position [m]')

subplot(212)
plot(t(1:10:end),ya_ddot_log(1:10:end),'b-',t(1:10:end),y_ddot_log(1:10:end),'r--', t(1:10:end),yddot_unpert(1:10:end),'k-.','LineWidth',2)
legend('ÿ_a','ÿ_c','ÿ_u')
axis([0 8 -20 20])
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')

