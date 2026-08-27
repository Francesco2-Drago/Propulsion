function [perimeter, area] = mesh_eval(Y, coordinate_system)

% INPUTS:
%   Y: mesh points in the specified coordinate system
%   coordinate_system: "cartesian" or "polar"
% OUTPUTS:
%   perimeter: perimeter of the mesh [m]
%   area: area of the mesh [m^2]

switch coordinate_system
    case "cartesian"
        x = Y(:,1);
        y = Y(:,2);
        perimeter = sum(hypot(diff(x), diff(y)));
        area = polyarea(x, y);          % area = abs(trapz(x, y));
    case "polar"
        theta = Y(:,1);
        rho = Y(:,2);
        dtheta = diff(theta);
        drho = diff(rho);
        perimeter = sum(hypot(drho, rho(1:end-1) .* dtheta));
        area = 0.5 * trapz(theta, rho.^2);
    otherwise
        error("Coordinate system must be 'cartesian' or 'polar'.");
end

% since we only created half geometry
perimeter = perimeter * 2;
area = area * 2;

end