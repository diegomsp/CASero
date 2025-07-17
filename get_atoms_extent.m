function extent = get_atoms_extent(atoms, edge_space)
% Returns the extent of the atoms along the x and y axes
%
% Args:
%   atoms: matrix of atomic positions (n_atoms x 3)
%   edge_space: amount of space to add on each side of the extent. Defaults to 0.

if nargin < 2
    edge_space = 0;
end

xmin = min(atoms(:, 1)) - edge_space;
xmax = max(atoms(:, 1)) + edge_space;
ymin = min(atoms(:, 2)) - edge_space;
ymax = max(atoms(:, 2)) + edge_space;

extent = [xmin, xmax, ymin, ymax];
end
