function main()
    % Case a
    alpha_A = 0.25; f_A = 0; beta_A = 0;
    [x_A1, M_A1] = RK4(alpha_A, f_A, beta_A, 0.7, 'A');
    [x_A2, M_A2] = RK4(alpha_A, f_A, beta_A, 1.5, 'A');

    % Case b
    alpha_B = 0; f_B = 0.005; beta_B = 0;
    [x_B1, M_B1] = RK4(alpha_B, f_B, beta_B, 0.7, 'B');
    [x_B2, M_B2] = RK4(alpha_B, f_B, beta_B, 1.5, 'B');

    % Case c
    alpha_C = 0; f_C = 0; beta_C = 50;
    [x_C1, M_C1] = RK4(alpha_C, f_C, beta_C, 0.5, 'C');
    [x_C2, M_C2] = RK4(alpha_C, f_C, beta_C, 2.0, 'C');

    % Plot Cases
    plotResults(x_A1, M_A1, x_A2, M_A2, 'Case A: Area Change');
    plotResults(x_B1, M_B1, x_B2, M_B2, 'Case B: Friction');
    plotResults(x_C1, M_C1, x_C2, M_C2, 'Case C: Heat Transfer');
 end

% RK4 implementation
function [x, M] = RK4(alpha, f, beta, M_0, caseType)
    x = linspace(0, 5, 21);
    h = x(2) - x(1);
    M = zeros(size(x));
    M(1) = M_0;

    for i = 1:length(x)-1
        if caseType == 'A'
            dMdx = caseA_ODE(x(i), M(i), alpha, f, beta);
        elseif caseType == 'B'
            dMdx = caseB_ODE(x(i), M(i), alpha, f, beta);
        elseif caseType == 'C'
            dMdx = caseC_ODE(x(i), M(i), alpha, f, beta);
        end

        k1 = h * dMdx;
        k2 = h * caseSpecificODE(x(i)+h/2, M(i)+k1/2, alpha, f, beta, caseType);
        k3 = h * caseSpecificODE(x(i)+h/2, M(i)+k2/2, alpha, f, beta, caseType);
        k4 = h * caseSpecificODE(x(i)+h, M(i)+k3, alpha, f, beta, caseType);

        M(i+1) = M(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

function dMdx= caseSpecificODE(x, M, alpha, f, beta, caseType)
    if caseType == 'A'
        dMdx = caseA_ODE(x, M, alpha, f, beta);
    elseif caseType == 'B'
        dMdx = caseB_ODE(x, M, alpha, f, beta);
    elseif caseType == 'C'
        dMdx = caseC_ODE(x, M, alpha, f, beta);
    end
end

function dMdx = caseA_ODE(x, M, alpha, f, beta)
    D = 1.0 + alpha*x;
    A = pi*D^2 / 4;
    term = - (pi * alpha * D) / (2 * A)
    dMdx = (M*(1+ 0.2*M^2)) / (1-M^2) * term;
end

function dMdx = caseB_ODE(x, M, alpha, f, beta)
    term = 0.5 * 1.4 * M^2 * (4*f/1.0);
    dMdx = (M * (1+0.2*M^2)) / (1-M^2) * term;
end

function dMdx = caseC_ODE(x, M, alpha, f, beta)
    T = 1000 + beta*x / 1.0;
    term = (1 + 1.4*M^2) * beta / (2* T * 1.0);
    dMdx = (M*(1 + 0.2*M^2)) / (1-M^2) * term;
end

% Plot function
function plotResults(x1, M1, x2, M2, caseName)
    figure;
    plot(x1, M1, 'b-o', x2, M2, 'r-s', 'LineWidth', 1.5);
    xlabel('Distance x (cm)');
    ylabel('Mach Number M');
    title(sprintf('%s\nM_0 = %.1f (blue) vs M_0 = %.1f (red)', caseName, M1(1), M2(1)));
    grid on;
    legend(sprintf('M_0=%.1f', M1(1)), sprintf('M_0=%.1f', M2(1)));
end
