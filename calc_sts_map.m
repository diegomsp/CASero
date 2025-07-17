function [final_map, extent] = calc_sts_map(cgeom, evecs, h, edge_space, dx, z_eff)

    evec = evecs;
    [orb_map, extent] = calc_orb_map(cgeom, evec, h, edge_space, dx, z_eff);
    final_map = abs(orb_map).^2;

end
