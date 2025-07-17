function [indices, local_grid] = get_local_grid(x_arr, y_arr, p, cutoff)
% Method that selects a local grid around an atom
%
% Args:
%   x_arr: global x array
%   y_arr: global y array
%   p: atomic position
%   cutoff (float, optional): extent of local grid in all directions. Defaults to 5.0.

if nargin < 4
    cutoff = 10.0;
end

x_min_i = find(abs(x_arr - p(1) + cutoff) == min(abs(x_arr - p(1) + cutoff)), 1);
x_max_i = find(abs(x_arr - p(1) - cutoff) == min(abs(x_arr - p(1) - cutoff)), 1);
y_min_i = find(abs(y_arr - p(2) + cutoff) == min(abs(y_arr - p(2) + cutoff)), 1);
y_max_i = find(abs(y_arr - p(2) - cutoff) == min(abs(y_arr - p(2) - cutoff)), 1);

local_x = repmat(x_arr(x_min_i:x_max_i), [y_max_i - y_min_i + 1, 1]);
local_y = repmat(y_arr(y_min_i:y_max_i)', [1, x_max_i - x_min_i + 1]);

indices = [x_min_i, x_max_i, y_min_i, y_max_i];
local_grid = {local_x, local_y};