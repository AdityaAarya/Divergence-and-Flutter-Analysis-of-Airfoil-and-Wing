
n = 2;                  % example structural or flow parameter
M = 0.6;                % Mach number (assumed constant)
Q = linspace(0.1, 5, 500);   % Varying Q

gamma_deg = [5, 10, 20, 30, 45];
colors = lines(length(gamma_deg));

figure;
hold on;
for i = 1:length(gamma_deg)
    gamma = deg2rad(gamma_deg(i));
    eta = zeros(size(Q));
    
    for j = 1:length(Q)
        q = Q(j);
        cosg = cos(gamma);
        sing = sin(gamma);
        root_term = sqrt(1 - M^2 * cosg^2);
        
        num = n * (1 - M^2 * cosg^2 - 2 * q * sing * root_term);
        den1 = q^2 * sing^2;
        den2 = (4 * n + 1) * q * sing * root_term;
        den3 = n * (1 - M^2 * cosg^2);
        
        eta(j) = num / (den1 + den2 + den3);
    end
    plot(Q, eta, 'LineWidth', 2, 'DisplayName', ['\gamma = ' num2str(gamma_deg(i)) '°'], ...
         'Color', colors(i, :));
end

xlabel('Q');
ylabel('\eta');
title('\eta vs Q for different sweep angles');
legend show;
grid on;
