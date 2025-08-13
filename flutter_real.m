clc;
clear;

% === Wing Parameters (Goland Wing) ===
l = 6.1;                % span (m)
b = 0.915;              % semi-chord (m)
EI = 9.765e6;           % bending stiffness (Nm²)
GJ = 0.989e6;           % torsional stiffness (Nm²)
m = 35.695;             % mass per unit length (kg/m)
I = 8.694;              % mass moment of inertia per unit length (kg·m)
rho = 1.225;            % air density (kg/m³)
a = -0.5;               % cg offset (dimensionless)

% === Sweep Angle (deg) ===
Lambda_deg = 30;                 % change this for different sweep angles
Lambda_rad = deg2rad(Lambda_deg);

% === Theodorsen Function Approximation ===
theodorsen_Ck = @(k) (1 - 0.165*1i*k) ./ (1 - 0.335*1i*k + 0.009*k.^2);

% === Generalized Mass and Stiffness Matrices ===
Mh = 0.5 * m * l;
Mt = 0.5 * I * l;
Mc = 0.5 * m * l * b * a;
M = [Mh, Mc; Mc, Mt];

Kh = (pi^4 * EI) / (2 * l^3);
Kt = 0.5 * GJ * (pi / l)^2;
K = [Kh, 0; 0, Kt];

% === Flutter Analysis ===
U_vals = linspace(5, 600, 1000);
max_real1 = zeros(size(U_vals));   % Real part of mode 1
max_real2 = zeros(size(U_vals));   % Real part of mode 2
flutter_U = NaN;

for idx = 1:length(U_vals)
    U = U_vals(idx);
    U_eff = U * cos(Lambda_rad);
    
    eigs_undamped = eig(M \ K);
    omega = sqrt(max(real(eigs_undamped)));
    k = omega * b / U_eff;
    Ck = theodorsen_Ck(k);
    
    % Aero damping and stiffness
    C_aero = 2 * pi * rho * b^2 * [0, U_eff; 0, (0.5 - a) * U_eff] * Ck;
    K_aero = 2 * pi * rho * b^2 * [0, U_eff^2; 0, (a + 0.5) * U_eff^2] * Ck;
    
    C_total = C_aero;
    K_total = K - K_aero;
    
    % State-space formulation
    A_upper = [zeros(2), eye(2)];
    A_lower = [-M \ K_total, -M \ C_total];
    A = [A_upper; A_lower];
    
    eigenvals = eig(A);
    real_parts = real(eigenvals);
    sorted_real = sort(real_parts, 'descend');
    
    max_real1(idx) = sorted_real(1);
    max_real2(idx) = sorted_real(2);
    
    if isnan(flutter_U) && any(real_parts >= 0)
        flutter_U = U;
    end
end

% === Plot Results ===
figure;
plot(U_vals, max_real1, 'b-', 'LineWidth', 1.8); hold on;
plot(U_vals, max_real2, 'r--', 'LineWidth', 1.8);
yline(0, 'k:');

if ~isnan(flutter_U)
    xline(flutter_U, '--k', sprintf('Flutter @ %.1f m/s', flutter_U), ...
          'LabelHorizontalAlignment', 'left', ...
          'LabelVerticalAlignment', 'bottom');
end

xlabel('Airspeed U (m/s)');
ylabel('Real Part of Eigenvalues');
title(['Flutter Analysis: Bending & Torsion Modes (\Lambda = ', num2str(Lambda_deg), '°)']);
legend('Mode 1 (Real Part)', 'Mode 2 (Real Part)', 'Location', 'northwest');
grid on;
