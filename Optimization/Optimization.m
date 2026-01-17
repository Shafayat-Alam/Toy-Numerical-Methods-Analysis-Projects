clear; clc;

x1 = 5000; %initial guess of x1               
alpha = 0.01;             
tol = 1e-5;               
max_iter = 1000;          

x1_history = x1;
z_history = compute_z(x1);

%steepest method iteration
for iter = 1:max_iter
    grad = compute_dzdx1(x1);     
    x1_new = x1 - alpha * grad;   

    x1_history(end+1) = x1_new;
    z_history(end+1) = compute_z(x1_new);

    if abs(x1_new - x1) < tol
        break;
    end

    x1 = x1_new;  
end

%Output of results
fprintf('x1: %.4f\n', x1);
x2 = (5092 - 1.040 * x1) / (0.0146 * x1 - 14);
x3 = (x1 + 19585.253) / 7.68;
fprintf('x2: %.4f\n', x2);
fprintf('x3: %.4f\n', x3);
fprintf('z: %.4f\n', compute_z(x1));

%z_value computation
function val = compute_z(x1)
    x2 = (5092 - 1.040 * x1) / (0.0146 * x1 - 14);
    x3 = (x1 + 19585.253) / 7.68;
    val = x1 + x2 + x3;
end

%steepest method computation
function grad = compute_dzdx1(x1)
    h = 1e-5;
    grad = (compute_z(x1 + h) - compute_z(x1 - h)) / (2 * h);
end
