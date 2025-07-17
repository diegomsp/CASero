function[filexyz] = rotate_to_XY(filexyz)

    u = G(2,:) - G(1,:);
    v = G(3,:) - G(1,:);
    uz=[0,0,1];
    e = cross( cross(u,v)/norm(cross(u,v)), uz); % axis of rotation
    alpha = atan2(norm(cross(e,uz)), dot(e,uz)); % angle of rotation

    for ii = 1:length(G)
        % Rodrigues rotation formula
        G(ii,:) = G(ii,:)*cos(alpha) + cross(e,G(ii,:))*sin(alpha) + dot(e,G(ii,:))*(1-cos(alpha))*e;
    end % ii

end % end function