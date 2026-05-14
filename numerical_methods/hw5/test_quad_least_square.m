% test_least_square.m

% --- Function 1: Good Approximation ---
f1 = @(x) exp(x);
x1 = linspace(-1, 1, 20)'; % Generate 20 points
y1 = f1(x1) + 0.1 * randn(size(x1)); % Add a tiny bit of noise for realism

[a1, b1, c1] = quad_least_square(x1, y1);
y1_approx = a1*x1.^2 + b1*x1 + c1;

% --- Function 2: Poor Approximation ---
f2 = @(x) sin(2*pi*x);
x2 = linspace(0, 2, 50)'; % Generate 50 points
y2 = f2(x2) + 0.1 * randn(size(x2)); % Add a tiny bit of noise

[a2, b2, c2] = quad_least_square(x2, y2);
y2_approx = a2*x2.^2 + b2*x2 + c2;

figure('Position', [100, 100, 1000, 400]);

% Subplot 1: Good Fit
subplot(1, 2, 1);
plot(x1, y1, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Data Points (e^x + noise)');
hold on;
plot(x1, y1_approx, 'b-', 'LineWidth', 2, 'DisplayName', 'Quadratic Fit');
title('Good Approximation (Exponential)');
xlabel('x'); ylabel('y');
legend('Location', 'northwest');
grid on;

% Subplot 2: Poor Fit
subplot(1, 2, 2);
plot(x2, y2, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Data Points (\sin(2\pi x) + noise)');
hold on;
plot(x2, y2_approx, 'r-', 'LineWidth', 2, 'DisplayName', 'Quadratic Fit');
title('Poor Approximation (Sine Wave)');
xlabel('x'); ylabel('y');
legend('Location', 'northeast');
grid on;

print('least_squares_comparison.png', '-dpng', '-r300');
disp('Plot saved as least_squares_comparison.png');