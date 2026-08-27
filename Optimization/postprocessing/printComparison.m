function printComparison(res)
    % printComparison
    % The comparative table across oxidizers, failures included.
    % INPUT
    %   res / struct / 1xM   Result records
    % OUTPUT
    %   None (prints to stdout)

    fprintf("\n%s\n", repmat('=', 1, 108));
    fprintf(" OXIDIZER COMPARISON\n");
    fprintf("%s\n", repmat('-', 1, 108));
    fprintf(" %-8s %-9s %9s %8s %7s %8s %7s %6s %8s %8s %7s  %s\n", ...
        "oxidizer", "geometry", "Isp_load", "mdot_ox", "O/F", "2R_c[m]", ...
        "L[m]", "L/D", "Gox0", "Gox_end", "sigma", "gate");
    fprintf("%s\n", repmat('-', 1, 108));
    for i = 1:numel(res)
        r = res(i);
        if r.ok
            fprintf(" %-8s %-9s %9.2f %8.2f %7.3f %8.3f %7.3f %6.2f %8.0f %8.1f %7.4f  %s\n", ...
                r.oxidizer, r.geometry, r.Isp_load, r.mdot_ox, r.OF_med, ...
                2*r.R_c, r.L, r.LD, r.Gox0, r.Gox_end, r.sigma, r.failTag);
        else
            fprintf(" %-8s %-9s %9s %8s %7s %8s %7s %6s %8s %8s %7s  %s\n", ...
                r.oxidizer, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", ...
                r.failTag);
        end
    end
    fprintf("%s\n", repmat('-', 1, 108));
end

