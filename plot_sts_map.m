function plot_sts_map(cgeom, evec, neighbor_list, filename, varargin)
    % Parámetros por defecto
    plot_bond = false;
    h = 10.0;
    edge_space = 5.0;
    dx = 0.1;
    z_eff = 3.25;

    % Verificar los argumentos de entrada opcionales
    if ~isempty(varargin)
        params = struct(varargin{:});
        if isfield(params, 'plot_bond')
            plot_bond = params.plot_bond;
        end
        if isfield(params, 'h')
            h = params.h;
        end
        if isfield(params, 'edge_space')
            edge_space = params.edge_space;
        end
        if isfield(params, 'dx')
            dx = params.dx;
        end
        if isfield(params, 'z_eff')
            z_eff = params.z_eff;
        end
    end

    % Obtener el mapa STS y su extensión
    [final_map, extent] = calc_sts_map(cgeom, evec, h, edge_space, dx, z_eff);

    % Graficar el mapa STS
    imagesc(extent(1:2), extent(3:4), final_map');
    colormap('hot');
    set(gca, 'YDir', 'normal');
    axis('image');

    % Graficar los enlaces atómicos si se especifica
    if plot_bond
        visualize_backbone(cgeom, neighbor_list);
    end
    caxis([0, max(max(final_map))*0.5]); % Ajustar brillo
    % Ajustar la figura al tamaño de la imagen y guardar
    set(gcf, 'PaperPositionMode', 'auto');
    print(filename, '-dpng', '-r0');
end



