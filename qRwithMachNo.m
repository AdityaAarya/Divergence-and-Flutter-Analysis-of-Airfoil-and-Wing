Kt = 1.2;
Clbeta0 = 0.05;
Cldelta0 = 0.08;
Cmacbeta0 = 0.03;
Cmacdelta0 = 0.04;
Clalpha0 = 5.7;
beta = 2*pi/180;
delta = 3*pi/180;
s = 20;
e1 = 0.9;
e2 = 1.1;
e3 = 0.7;
c = 3;
M = linspace(0.1, 0.7, 50);
qR = zeros(size(M));
for i = 1:length(M)
    mach = M(i);
    numerator = Kt * (Clbeta0 * beta + Cldelta0 * delta);
    denominator1 = s * Clbeta0 * (e1 * (Clalpha0 / sqrt(1 - mach^2)) + e2)*beta;
    denominator2 = s * Cldelta0 * (e1 * (Clalpha0 / sqrt(1 - mach^2)) - e3)*delta;
    denominator3 = s * c * (Cmacbeta0 * beta + Cmacdelta0 * delta);
    qR(i) = numerator / (denominator1 + denominator2 - denominator3);
end

figure;
plot(M, qR, 'b', 'LineWidth', 2);
xlabel('Mach number (M)');
ylabel('q_R');
title('Variation of q_R with Mach number');
grid on;
