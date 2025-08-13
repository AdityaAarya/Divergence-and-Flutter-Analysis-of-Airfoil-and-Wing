clc;
clear;

% flutter whatever
c=0.2540;
b = c/2;              % semi-chord
m = 3.3843;             % mass per unit length
I = 0.0135;              % mass moment of inertia
rho = 1.225;            % density
a = -0.5;               % offset
sa=0.0859;              %instability
Kh = 2818.8;               % Bending stiffness (N/m)
Kt = 37.3;               % Torsional stiffness (Nm/rad)

% sweep if needed
Lambda_deg = 0;
Lambda_rad = deg2rad(Lambda_deg);

% jones approximation valid in all incompressibles
theodorsen_Ck = @(k) (1 - 0.165*1i*k) ./ (1 - 0.335*1i*k + 0.009*k.^2);
% change in case of continuous EI or GJ whatever
Mh = m ;
Mt = I;
Mc = sa;
M = [Mh, Mc; Mc, Mt];
K = [Kh, 0; 0, Kt];

U_vals = linspace(1, 30, 1000);
max_real1 = zeros(size(U_vals));
max_real2 = zeros(size(U_vals));
flutter_U = NaN;

for idx = 1:length(U_vals)
    U = U_vals(idx);
    U_eff = U * cos(Lambda_rad);
    M_nc=pi*rho*b^2*[1,-b*a;-b*a,(0.125+a^2)*b^2];
    eigs_undamped = eig(M+(M_nc/U_eff) \ K);
    omega = sqrt(max(real(eigs_undamped)));
    k = omega * b / U_eff;
    Ck = theodorsen_Ck(k);
    
    % Aero damping and stiffness non circulatory and circulatory
    C_nc=[0,pi*rho*U_eff*b^2;0,pi*rho*U_eff*(0.5-a)*b^3];
    M_nc=pi*rho*b^2*[1,-b*a;-b*a,(0.125+a^2)*b^2];
    C_aero = pi * rho * b^2 * [2/b, 2*U_eff*(0.5-a); -2*U_eff*(a+0.5), -2*U_eff*(a+0.5)*(0.5-1)*b] * Ck;
    K_aero = 2 * pi * rho * b * [0, U_eff^2; 0, -(a + 0.5)*b* U_eff^2] * Ck;
    M_total=M+M_nc;
    C_total = C_aero+C_nc;
    K_total = K + K_aero;
    
    % EVP
    A_upper = [zeros(2), eye(2)];
    A_lower = [-M_total \ K_total, -M_total \ C_total];
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
title(['Flutter Analysis: Input Kh & Kt, Sweep \Lambda = ', num2str(Lambda_deg), '°']);
legend('Mode 1 (Real Part)', 'Mode 2 (Real Part)', 'Location', 'northwest');
grid on;
