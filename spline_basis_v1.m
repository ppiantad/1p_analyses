%% Set up parameters for the spline basis
n_basis = 25;          % Number of basis functions
n_order = 4;           % Order of the spline (degree + 1)

% Create B-spline basis with 200 rows
n_rows = 200;  % Adjust as needed
basis_set = zeros(n_rows, n_basis);

% Create and evaluate basis functions
for i = 1:n_rows
    % Create original break points with 27 inner points
    original_break_points = linspace(0, 1, 27);
    
    % Adjust break points to fit the window size for the current row
    break_points_adjusted = linspace(0, 1, 27);  % Ensure original break points are used
    break_points_adjusted = (i - 1) / (n_rows - 1) + break_points_adjusted / (n_rows - 1);
    
    % Make sure the highest value is exactly 1
    break_points_adjusted(end) = 1;
    
    % Create B-spline basis
    basis = create_bspline_basis([0, 1], 27, n_order, break_points_adjusted);
    
    % Evaluate basis functions at specific time points
    eval_points = linspace(0, 1, 200);  % Adjust as needed
    basis_functions = eval_basis(eval_points, basis);
    
    % Extract the middle 25 columns
    basis_set(i, :) = basis_functions(:, 2:26);
end