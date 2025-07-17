function [orb_map, extent] = calc_orb_map(cgeom, evec, h, edge_space, dx, z_eff)
% Calculates the molecular orbital map using the given geometry and eigenvectors
%
% Args:
%   cgeom: matrix of atomic positions (n_atoms x 3)
%   evec: vector of eigenvectors (n_atoms x 1)
%   h (float, optional): extent of the orbital in the z direction. Defaults to 10.0.
%   edge_space (float, optional): amount of space to add on each side of the extent. Defaults to 5.0.
%   dx (float, optional): grid spacing. Defaults to 0.1.
%   z_eff (float, optional): effective nuclear charge interacting with the pz orbital. Defaults to 3.25.

%%%%%%%%%%Nota: en la función get_local_grid, el output es una lista [ [x_min_i, x_max_i, y_min_i, y_max_i],
% [local_x, local_y] ], y se utiliza la sintaxis de indexación de lista para extraer los valores de los 
% índices [0], [1], [2], [3] en Python. En Matlab, el output es una tupla 
% ( [x_min_i, x_max_i, y_min_i, y_max_i], {local_x, local_y} ), y se 
% utiliza la sintaxis de indexación de tupla para extraer los valores [1]{1}, [1]{2} en lugar de [1][0], [1][1].

%if nargin < 6
%    z_eff = 3.25;
%end
%if nargin < 5
%    dx = 0.1;
%end
%if nargin < 4
%    edge_space = 5.0;
%end
%if nargin < 3
%    h = 10.0;
%end

extent = get_atoms_extent(cgeom, edge_space);

% define grid
x_arr = extent(1):dx:extent(2);
y_arr = extent(3):dx:extent(4);

% update extent so that it matches with grid size
extent(2) = x_arr(end) + dx;
extent(4) = y_arr(end) + dx;

orb_map = zeros(length(x_arr), length(y_arr));

for i = 1:size(cgeom, 1)
    at = cgeom(i, :);
    coef = evec(i);
    p = at;
    [local_i, local_grid] = get_local_grid(x_arr, y_arr, p, 1.2*h + 4.0);
    pz_orb = carbon_2pz_slater(local_grid{1} - p(1), local_grid{2} - p(2), h, z_eff);
    orb_map(local_i(1):local_i(2), local_i(3):local_i(4)) = orb_map(local_i(1):local_i(2), local_i(3):local_i(4)) + coef*pz_orb';
end

end
