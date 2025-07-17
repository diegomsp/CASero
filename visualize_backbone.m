function visualize_backbone(atoms, neighbor_list)
    i_arr = neighbor_list(1, :);
    j_arr = neighbor_list(2, :);
    for k = 1:numel(i_arr)
        i = i_arr(k);
        j = j_arr(k);
        if i < j
            p1 = atoms(i, :);
            p2 = atoms(j, :);
            plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', 3.0, 'LineJoin', 'round');
        end
    end
end

