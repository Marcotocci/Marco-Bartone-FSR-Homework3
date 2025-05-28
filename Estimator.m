clear 
close all
clc

%% Quadrotor data

data = load('ws_homework_3_2025.mat');

%Parameters

m = 1.5;

Ib = diag([1.2416, 1.2416, 2*1.2416]);

g = 9.81;

e_3 = [0; 0; 1];

%Data
t = data.linear_vel.time;

u_T = data.thrust.signals.values;

tau_b = data.tau.signals.values;

pb_dot = data.linear_vel.signals.values;

eta = data.attitude.signals.values;

eta_dot = data.attitude_vel.signals.values;

Ts = 0.001;

r = 1; %estimator order

q = zeros(6,length(t));

gamma = zeros(6,length(t),r);

ext_wrench_e = zeros(6,length(t));

%% Estimator

k0 = 100;

s = tf('s');

G = k0^r/(s + k0)^r;

c = cell2mat(G.Denominator);
c = c(2:end);

k = zeros(r,1);

x = zeros(6,length(t));
y = zeros(6,length(t));
x_next = zeros(6,length(t),r);
y_next = zeros(6,length(t),r);

% Computing ki
P = 1;
for j = 1:r
    k(j) = c(j)/P;
    P = P*k(j);
end
k = flip(k);

% Loop 
for i = 1:size(t,1)-1
    % computing Rb, Q, Q_dot, C at time i
    Rb = attitude(eta(i,1), eta(i,2), eta(i,3));

    [Q, Q_dot] = Q_generator(eta(i,1), eta(i,2), eta_dot(i,1), eta_dot(i,2));

    C = Centrif_Matrix(eta_dot(i,:), Q, Q_dot, Ib);

    %compute q

    M = Q'*Ib*Q;

    q(:,i+1) = [m*eye(3), zeros(3,3); zeros(3,3), M]*[pb_dot(i+1,:)'; eta_dot(i+1,:)'];

    x(:,i) = ext_wrench_e(:,i) + [(m*g*e_3 - u_T(i)*Rb*e_3); C'*eta_dot(i,:)' + Q'*tau_b(i,:)'];

    y(:,i+1) = y(:,i) +Ts*x(:,i);

    %gamma

    gamma(:,i+1,1) = k(1) * (q(:,i+1) - y(:,i+1));

    % loop for gamma
   for j=2:r
       x_next(:,i,j) = -ext_wrench_e(:,i) + gamma(:,i,j-1);
       y_next(:,i+1,j) = y_next(:,i,j) + Ts* x_next(:,i,j);
       gamma(:,i+1,j) = k(j) * y_next(:,i+1,j);
   end

   % update estimation

   ext_wrench_e(:,i+1) = gamma(:,i+1,r);

end

%% Plot

y_label_forces = {'$f_x$ [N]', '$f_y$ [N]', '$f_z$ [N]'};
y_label_torques = {'$\tau_x$ [Nm]', '$\tau_y$ [Nm]', '$\tau_z$ [Nm]'};
x_label = '$t$ [s]';

%Expected values

f_x = 1; 
f_y = 1;
tau_z = -0.4;

fig_forces = latex_triple_subplot_plot(t, ...
    ext_wrench_e(1,:), ext_wrench_e(2,:), ext_wrench_e(3,:), ...
    y_label_forces{1}, y_label_forces{2}, y_label_forces{3}, ...
    x_label, ...
    'Estmation external forces');

subplot(3,1,1); hold on;
plot(t, f_x*ones(size(t)), 'r--', 'LineWidth', 2);
legend({'$\hat{f}_x$', '$f_x$'}, 'Interpreter', 'latex', 'Location', 'southeast');

subplot(3,1,2); hold on;
plot(t, f_y*ones(size(t)), 'r--', 'LineWidth', 2);
legend({'$\hat{f}_y$', '$f_y$'}, 'Interpreter', 'latex', 'Location','southeast');

exportgraphics(fig_forces,'Estimator_forces_r1_k100.pdf');


fig_torques = latex_triple_subplot_plot(t, ...
    ext_wrench_e(4,:), ext_wrench_e(5,:), ext_wrench_e(6,:), ...
    y_label_torques{1}, y_label_torques{2}, y_label_torques{3}, ...
    x_label, ...
    'Estmation external torques');

subplot(3,1,3); hold on;
plot(t, tau_z*ones(size(t)), 'r--', 'LineWidth', 2);
legend({'$\hat{\tau}_z$', '$\tau_z$'}, 'Interpreter', 'latex');

exportgraphics(fig_torques, 'Estimator_torques_r1_k100.pdf');


%% Compute mass of the UAV from estimated disturbances

m_tilde = ext_wrench_e(3,end)/g;
m_real = m + m_tilde;


%% Functions

function Rb = attitude(phi, theta, psi)

    Rb = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(psi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
          cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
              -sin(theta)    ,                     sin(phi)*cos(theta),                               cos(phi)*cos(theta)       ];

end

function [Q, Q_dot] = Q_generator(phi, theta, phi_d, theta_d)

    Q = [1, 0, -sin(theta);
        0, cos(phi), cos(theta)*sin(phi);
        0, -sin(phi), cos(theta)*cos(phi)];

    Q_dot = [0, 0, -cos(theta)*theta_d;
             0, -sin(phi)*phi_d, -sin(theta)*theta_d*sin(phi) + cos(theta)*cos(phi)*phi_d;
             0, -cos(phi)*phi_d, -sin(theta)*theta_d*cos(phi) - cos(theta)*sin(phi)*phi_d];

   
end

function S =skew(w)

    S = [0, -w(3), w(2);
         w(3), 0, -w(1);
         -w(2), w(1), 0];
end

function C = Centrif_Matrix(eta_dot, Q, Q_dot, Ib)
    
    S = skew(Q*eta_dot');

    C = Q'*S*Ib*Q + Q'*Ib*Q_dot;

end
