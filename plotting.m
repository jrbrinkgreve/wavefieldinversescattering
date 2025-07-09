% Parameters
Nx = 21;
Ny = 21;
x = 1:Nx;
y = 1:Ny;
[X, Y] = meshgrid(x, y);
source_location = [11, 201];  % in lambda/2 units

% Receiver locations (in original coord system)
M_per_array = 20;
rec_locations_x_original = linspace(-19, 41, M_per_array);
rec_locations_y_original = 31 * ones(1, M_per_array);

% Plot
figure;
hold on;
axis equal;

% Imaging domain
scatter(X(:), Y(:), 30, 'b', 'filled');  % Imaging grid

% Source
plot(source_location(1), source_location(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
%text(source_location(1) + 1, source_location(2), 'Source', 'Color', 'r');

% Receivers
plot(rec_locations_x_original, rec_locations_y_original, 'g-', 'LineWidth', 2);
%text(mean(rec_locations_x_original), rec_locations_y_original(1) + 2, 'Receiver Array', 'HorizontalAlignment', 'center', 'Color', 'g');



% Labels and formatting
xlabel('x [\lambda / 20]');
ylabel('y [\lambda / 20]');
title('Imaging Configuration');
legend('Imaging Grid', 'Source', 'Receivers');
grid on;
xlim([-25 45]);
ylim([-20 210]);
set(gca, 'FontSize', 20);  % axes tick labels

% Make y-axis point down
set(gca, 'YDir', 'reverse');





%-----------------------

% Plot
figure;
hold on;
axis equal

% Imaging domain
scatter(X(:), Y(:), 30, 'b', 'filled');  % Imaging grid

% Source
plot(source_location(1), source_location(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Receiver arrays as green lines
plot(rec_locations_x_original, rec_locations_y_original, 'g-', 'LineWidth', 2); % top
plot(rec_locations_x_west, rec_locations_y_west, 'g-', 'LineWidth', 2);         % left
plot(rec_locations_x_east, rec_locations_y_east, 'g-', 'LineWidth', 2);         % right
plot(rec_locations_x_north, rec_locations_y_north, 'g-', 'LineWidth', 2);       % bottom

% Labels and formatting
xlabel('x [\lambda / 20]');
ylabel('y [\lambda / 20]');
title('Imaging Configuration');
legend('Imaging Grid', 'Source', 'Receivers');
grid on;
xlim([-25 45]);
ylim([-20 210]);
set(gca, 'FontSize', 20);  % axes tick labels

% Make y-axis point down
set(gca, 'YDir', 'reverse');



%{
% Labels and formatting
xlabel('x [\lambda / 20]', 'FontSize', 20);
ylabel('y [\lambda / 20]', 'FontSize', 20);
title('Sensor, Source, and Imaging Domain Configuration', 'FontSize', 22);
legend('Imaging Grid', 'Source', 'Receivers');
grid on;
xlim([-25 45]);
ylim([-40 210]);
set(gca, 'FontSize', 20);

% Make y-axis point down
set(gca, 'YDir', 'reverse');

%}
