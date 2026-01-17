clc; clear;

% ======= Stage 1: Solve Equations 1,2, and 6 to find CP, NHR, NKW =======

x = [1.6; 6800; 360000];% Initial Guess Vector: [CP, NHR, NKW]

% Convergence settings
tolerance = 1e-12; 
max_iter = 1000;
WFR = 145000; %WFR is given

for iter = 1:max_iter
    CP = x(1);
    NHR = x(2);
    NKW = x(3);

    HL = ((NHR - 3412)*NKW)/1e6; % HL from Equation 6

    % Forming Function Vector
    F = zeros(3,1);
    F(1) = -45.19*CP^4 + 420*CP^3 - 1442*CP^2 + 2248*CP + 6666 - NHR;
    F(2) = 4883*CP^4 - 44890*CP^3 + 152600*CP^2 - 231500*CP + 383400 - NKW;
    F(3) = ((NHR - 3412)*NKW)/1e6 - HL;  

    J = zeros(3,3); % --- Jacobian Matrix

    % Jacobi for Equation 1
    J(1,1) = -180.76*CP^3 + 1260*CP^2 - 2884*CP + 2248;
    J(1,2) = -1;
    J(1,3) = 0;

    % Jacobi for Equation 2
    J(2,1) = 19532*CP^3 - 134670*CP^2 + 305200*CP - 231500;
    J(2,2) = 0;
    J(2,3) = -1;

    % Jacobi for Equation 3
    J(3,1) = 0;
    J(3,2) = NKW / 1e6;
    J(3,3) = (NHR - 3412) / 1e6;

    % Newton-Raphson 
    dx = -J \ F;
    x_new = x + dx;

    % Checking Relative error
    rel_error = abs((x_new - x) ./ x_new) * 100;

    % Display status
    %fprintf('Iter %3d | CP = %.4f | NHR = %.2f | NKW = %.2f | Max RelErr = %.4e\n', ...
     %   iter, x(1), x(2), x(3), max(rel_error));

    if all(rel_error < tolerance)
        disp('Stage 1 Converged.');
        break;
    end

    x = x_new; % Setting new x for iteraton
end

% Final Outputs
CP = x(1);
NHR = x(2);
NKW = x(3);
HL = ((NHR - 3412)*NKW)/1e6;

fprintf('\nFinal Stage 1 Results:\n');
fprintf('CP  = %.6f\n', CP);
fprintf('NHR = %.2f\n', NHR);
fprintf('NKW = %.2f\n', NKW);
fprintf('HL  = %.2f\n', HL);

% ===== Stage 2: Direct Calculation of CR using Equation 5 ======= 

% HL from Stage 1
% WFR, given Water Flow Rate
CR = (2000 * HL) / WFR;  % Equation 5

fprintf('\nFinal Stage 2 Results (Updated):\n');
fprintf('CR  = %.2f\n', CR);

% ======= Stage 3: Solve Equation 3 for CWT  =========

% CP and HL known from Stage 1;
% WFR is given 

CWT = 85; % Initial guess for CWT

% Convergence settings
tolerance = 1e-12;
max_iter = 1000;

for iter = 1:max_iter
    % Equation 3
    f = 1.6302 ...
        - 0.050095 * CWT ...
        + 0.00055796 * CWT^2 ...
        + 0.00032946 * HL ...
        - 0.000010229 * HL * CWT ...
        + 0.00000016253 * HL * CWT^2 ...
        + 0.00000042658 * HL^2 ...
        - 0.000000092331 * HL^2 * CWT ...
        + 0.000000000071265 * HL^2 * CWT^2 ...
        - CP;

    % Equation 3 Jacobian
    df = -0.050095 ...
         + 2 * 0.00055796 * CWT ...
         - 0.000010229 * HL ...
         + 2 * 0.00000016253 * HL * CWT ...
         - 0.000000092331 * HL^2 ...
         + 2 * 0.000000000071265 * HL^2 * CWT;

    % Newton-Raphson 
    dCWT = -f / df;
    CWT_new = CWT + dCWT;

    rel_err = abs((CWT_new - CWT) / CWT_new) * 100; % Relative error check

    if rel_err < tolerance
        disp('Stage 3 Converged.');
        break;
    end

    CWT = CWT_new; % Newton-Raphson iteration
end

% Final Output
fprintf('\nFinal Stage 3 Result:\n');
fprintf('CWT = %.2f\n', CWT);

% ===== Stage 4: Direct Calculation of WBT from Equation 4 ======

% Using known values CR, CWT, WFR 

% Rearrange Equation 4 to solve for WBT: CWT = A + B*WFR + C*CR + D*CR*WFR + E*WBT + F*WBT*WFR + G*WBT*CR + H*WBT*CR*WFR
%                                        WBT = (CWT - (A + B*WFR + C*CR + D*CR*WFR)) / (E + F*WFR + G*CR + H*CR*WFR)

% Coefficients of Equation 4
A = -10.046;
B = 0.00022801;
C = 0.85396;
D = 0.0000018617;
E = 1.0957;
F = -0.000022425;
G = -0.011978;
H = 0.00000014378;

% Calculate WBT 
numerator = CWT - (A + B*WFR + C*CR + D*CR*WFR);
denominator = E + F*WFR + G*CR + H*CR*WFR;
WBT = numerator / denominator;

% Final Output
fprintf('\nStage 4 Result:\n');
fprintf('WBT = %.2f\n', WBT);

AP = CWT - WBT;

% Results Summary
fprintf('\n========== Final Results Summary ==========\n');
fprintf('CP  = %.6f\n', CP);
fprintf('NHR = %.2f\n', NHR);
fprintf('NKW = %.2f\n', NKW);
fprintf('HL  = %.2f\n', HL);

fprintf('WFR = %.0f GPM (Given)\n', WFR);
fprintf('CR  = %.2f\n', CR);

fprintf('CWT = %.2f\n', CWT);

fprintf('WBT = %.2f\n', WBT);
fprintf('Tower Approach (AP) = %.2f\n', AP)
fprintf('===========================================\n');
