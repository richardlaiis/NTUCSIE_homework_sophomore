% test_spline.m

f = @(x) sin(x);
df = @(x) cos(x);

x = linspace(0, 2*pi, 6);
xq = linspace(0, 2*pi, 200);

y_vals = f(x);
y_with_slopes = [df(x(1)), y_vals, df(x(end))];

yq_builtin = spline(x, y_with_slopes, xq);
yq_custom  = my_spline(x, y_with_slopes, xq);

max_err = max(abs(yq_builtin - yq_custom));
fprintf('Maximum difference between custom and built-in: %e\n', max_err);

fig = figure('visible', 'off');

set(fig, 'Position', [0, 0, 800, 600]);

plot(x, y_vals, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Data Points');
hold on;
plot(xq, f(xq), 'k-', 'LineWidth', 1, 'DisplayName', 'True Function \sin(x)');
plot(xq, yq_builtin, 'b-', 'LineWidth', 4, 'DisplayName', 'Built-in spline');
plot(xq, yq_custom, 'r--', 'LineWidth', 2, 'DisplayName', 'my\_spline');

legend('Location', 'northeastoutside');
title('Cubic Spline: Slope Boundary Conditions');
xlabel('x');
ylabel('f(x)');
grid on;
hold off;

output_filename = 'spline_plot.png';
print(fig, output_filename, '-dpng', '-r300');
fprintf('Plot successfully saved to %s\n', output_filename);