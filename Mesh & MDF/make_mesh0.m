function Y = make_mesh0(geometry_type, meshdata, npoints, coordinate_system)

% INPUTS:
%   geometry_type: "star", "cylinder"
%   meshdata: struct with variables for geometry calculation
%   npoints: number of points in the mesh
%   coordinate_system: "cartesian" or "polar"
% OUTPUTS:
%   Y: mesh points in the specified coordinate system

switch geometry_type
    case "star"
        ri = meshdata.inner_diameter/2;
        re = meshdata.outer_diameter/2;
        n_tips = meshdata.n_tips;

        % ATTENZIONE ZIOPERA! QUESTO PASSAGGIO IL SANTOLINI NON LO FA
        % ri is the apothem of the inner polygon, so the polygon vertices are at ri/cos(pi/n_tips)
        r_vertex = ri / cos(pi/n_tips);
        % r_vertex = ri;

        % upper-boundary star vertices (theta from 0 to pi), alternating tip and polygon vertex
        theta_vertices = (0:n_tips)' * (pi/n_tips);
        rho_vertices = zeros(n_tips + 1, 1);
        rho_vertices(1:2:end) = re;
        rho_vertices(2:2:end) = r_vertex;

        [x_vertices, y_vertices] = pol2cart(theta_vertices, rho_vertices);

        % Sample points uniformly along edge length to keep all sides straight
        ds = hypot(diff(x_vertices), diff(y_vertices));
        s_vertices = [0; cumsum(ds)];
        s_query = linspace(0, s_vertices(end), npoints)';
        x = interp1(s_vertices, x_vertices, s_query, "linear");
        y = interp1(s_vertices, y_vertices, s_query, "linear");

    case "cylinder"
        diam = meshdata.diameter;
        theta = linspace(0, pi, npoints)';
        rho = diam/2 * ones(npoints, 1);
        [x, y] = pol2cart(theta, rho);
    
    case "rod_and_tube"
        rod_diam = meshdata.inner_diameter;
        tube_diam = meshdata.outer_diameter;

        theta_rod = linspace(0, pi, npoints)';
        theta_tube = linspace(pi, 0, npoints)';
        rho_rod = rod_diam/2 * ones(npoints, 1);
        rho_tube = tube_diam/2 * ones(npoints, 1);
        [x_rod, y_rod] = pol2cart(theta_rod, rho_rod);
        [x_tube, y_tube] = pol2cart(theta_tube, rho_tube);
        x = [x_rod; x_tube];
        y = [y_rod; y_tube];
    
    otherwise
        error("Geometry type must be 'star' or 'cylinder'.");
end

switch coordinate_system
    case "cartesian"
        Y = [x, y];
    case "polar"
        [theta, rho] = cart2pol(x, y);
        theta = mod(theta, 2*pi);
        Y = [theta, rho];
    otherwise
        error("Coordinate system must be 'cartesian' or 'polar'.");
end

end