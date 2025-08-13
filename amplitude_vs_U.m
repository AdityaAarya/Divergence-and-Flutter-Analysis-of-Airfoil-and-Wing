clc;
clear;

% Given constants
c = 0.2540;
b = c/2;
m = 3.3843;
I = 0.0135;
rho = 1.225;
a = -0.5;
sa = 0.0859;
Kh = 2818.8;
Kt = 37.3;

Mh = m;
Mt = I;
Mc = sa;
M = [Mh, Mc; Mc, Mt];
K = [Kh, 0; 0, Kt];

% Velocity values
U_vals = [5, 15, 30];  % representative velocities to simulate
tspan = [0 5];         % simulation time span
init_cond = [0.01; 0.01; 0; 0]; % [h; theta; h_dot; theta_dot] initial conditions
theodorsen_Ck = @(k) besselh(1,2,1i*k) ./ (besselh(1,2,1i*k) + 1i * besselh(0,2,1i*k));
figure;
for i = 1:length(U_vals)
    U_eff = U_vals(i);

    % Compute reduced frequency and Theodorsen coefficient
    eigs_undamped = eig(M + (1/U_eff)*(M \ K));
    omega = sqrt(max(real(eigs_undamped)));
    k = omega * b / U_eff;
    Ck = theodorsen_Ck(k);

    % Non-conservative Mass and Damping matrices
    M_nc = pi*rho*b^2 * [1, -b*a; -b*a, (0.125 + a^2)*b^2];
    C_nc = [0, pi*rho*U_eff*b^2; 
            0, pi*rho*U_eff*(0.5 - a)*b^3];

    % Aerodynamic damping and stiffness
    C_aero = pi * rho * b^2 * ...
        [2/b, 2*U_eff*(0.5 - a); 
        -2*U_eff*(a + 0.5), -2*U_eff*(a + 0.5)*(0.5 - 1)*b] * Ck;

    K_aero = 2 * pi * rho * b * ...
        [0, U_eff^2; 
         0, -(a + 0.5)*b*U_eff^2] * Ck;

    % Total matrices
    M_total = M + M_nc;
    C_total = C_nc + C_aero;
    K_total = K + K_aero;

    % Convert to state-space (1st-order system)
    A = [zeros(2), eye(2);
        -M_total \ K_total, -M_total \ C_total];

    % Solve ODE
    ode_sys = @(t, x) A * x;
    [t, x] = ode45(ode_sys, tspan, init_cond);

    % Plot results
    subplot(length(U_vals), 1, i);
    plot(t, x(:,1), 'b', 'LineWidth', 1.5); hold on;
    plot(t, x(:,2), 'r', 'LineWidth', 1.5);
    title(['Oscillatory Response at U = ', num2str(U_eff), ' m/s']);
    xlabel('Time (s)');
    ylabel('Displacement');
    legend('h(t)', '\theta(t)');
    grid on;
end
