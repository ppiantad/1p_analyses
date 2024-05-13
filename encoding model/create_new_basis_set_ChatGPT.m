% Import or load your data into MATLAB

% Step 1: Define the range for the basis functions
min_time = 0;  % Minimum time point
max_time = 21;  % Maximum time point

% Step 2: Initialize the basis object
nbasis = 21;  % Number of basis functions (splines)
ntpts = 200;  % Number of time points

% Create a basis object with nbasis B-spline basis functions over the
% interval [min_time, max_time]
basisobj = create_bspline_basis([min_time, max_time], nbasis+4, ntpts, 0:nbasis);

% Optionally, you can set other parameters such as the range of data, knots, etc.

% Example: specifying knots
% norder = 4; % order of B-spline
% basisobj = create_bspline_basis([min_time, max_time], nbasis, ntpts, norder);

% You can further customize the basis object as per your requirements

% Display the basis functions
plot(basisobj);

% Optionally, you can use the basis object for further analysis

% Define the time points at which you want to evaluate the basis functions
eval_time = linspace(min_time, max_time, ntpts);

% Get the basis matrix evaluated at the specified time points
basis_matrix = getbasismatrix(eval_time, basisobj);

% Convert the sparse basis matrix to a full double array
basis_matrix_full = full(basis_matrix);