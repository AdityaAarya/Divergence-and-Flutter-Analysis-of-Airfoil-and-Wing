% Parameters
c = 0.127;      % m
l = 0.762;      % m
e = 0.25 * c;   % m
EI = 25.3232;   % N-m^2
GJ = 38.22;     % N-m^2
rho = 1.225;    % kg/m^3
AR = (2 * l) / c;
Clalpha0 = 2 * pi;

% Sweep angles
gamma_deg = linspace(-45, 0, 100);
gamma_rad = deg2rad(gamma_deg);

% Mach number range
M_vals = linspace(0.1, 0.7, 7);

% Prepare figure
figure;
hold on;

% Loop over Mach numbers
for i = 1:length(M_vals)
    M = M_vals(i);

    % Compute Clalpha (compressibility and sweep corrected)
    beta = sqrt(1 - (M^2) * (cos(gamma_rad)).^2);
    term1 = AR * (Clalpha0 * cos(gamma_rad)) ./ beta;
    term2 = AR * sqrt(1 + ((Clalpha0 * cos(gamma_rad)).^2) ./ (AR * pi * beta).^2);
    Clalpha = term1 ./ term2 + (1/pi) * beta;

    % Compute qD and UD
    numerator = GJ * (pi^2);
    correction_term = 1 - ((0.03515 * GJ .* tan(gamma_rad) * l * (pi^2)) ./ (e * EI));
    denominator = (4 * (l^2) .* (cos(gamma_rad).^2) * e * c .* Clalpha) .* correction_term;

    qD = numerator ./ denominator;
    UD = sqrt((2 * qD) / rho);

    % Plot
    plot(gamma_deg, UD, 'DisplayName', ['M = ' num2str(M)], 'LineWidth', 1.5);
end

xlabel('Sweep Angle (degrees)');
ylabel('Divergence Speed (m/s)');
title('Divergence Speed vs Sweep Angle for Varying Mach Numbers');
grid on;
legend('Location', 'best');
hold off;
