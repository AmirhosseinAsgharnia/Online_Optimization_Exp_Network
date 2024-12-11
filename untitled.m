% Generate example data
% data = rand(1000, 2);
% 
% % Plot histogram
% histogram2(data, 'Normalization', 'pdf');
% hold on;

% % Perform Kernel Density Estimation
% [x_vals, kde_vals] = ksdensity(data,0:0.1:1);
% 
% % Plot KDE
% plot(kde_vals,x_vals, 'LineWidth', 2, 'Color', 'r');
% legend('Histogram', 'KDE');
% hold off;

% Generate example data
x = randn(1000, 1);
y = randn(1000, 1);

% Define number of bins
num_bins = 30;

% Create 2D histogram
[counts, x_edges, y_edges] = histcounts2(x, y, num_bins, num_bins);

% Create grid for plotting
[x_grid, y_grid] = meshgrid(x_edges(1:end-1), y_edges(1:end-1));

% Display using surf
surf(x_grid, y_grid, counts', 'EdgeColor', 'none');
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Counts');
title('2D Histogram - Surface Plot');