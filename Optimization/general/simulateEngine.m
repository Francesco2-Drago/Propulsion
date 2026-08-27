function r = simulateEngine(R_c, L, At, mdot_ox, ctx, lookupN, cfg)
    % simulateEngine
    % Run the exact ODE combustion model on a normalized geometry lookup
    % rescaled to a given casing radius. This is the single evaluator of a
    % candidate engine: sizeEngine calls it to solve C6, C7 and C8, and
    % verifySolver calls it to sweep one knob at a time, so the solver and its
    % verification cannot drift apart.
    %
    % The burn is stopped by whichever comes first: the model's own events
    % inside run_combustion (O/F out of the CEA domain, pressure, burnout) or
    % the C12 residual web, which is imposed here through b_max.
    % INPUT
    %   R_c     / double / 1x1   Casing radius [m]
    %   L       / double / 1x1   Grain length [m]
    %   At      / double / 1x1   Throat area [m2]
    %   mdot_ox / double / 1x1   Oxidizer mass flow [kg/s]
    %   ctx     / struct / 1x1   Context from engineContext (thermo, engine, K)
    %   lookupN / struct / 1x1   Normalized geometry lookup at R_c = 1
    %   cfg     / struct / 1x1   n, a, rho_f, web_min [m], x2 [-], p_min, p_max
    %                            [Pa], time_output [s], fine_ode
    % OUTPUT
    %   r       / struct / 1x1   run_combustion output; r.ok false with r.err on
    %                            failure

    K = ctx.K;

    % Rescale the normalized lookup. Areas go as R_c^2, lengths as R_c, and
    % dperim_db is scale invariant.
    lk.b = R_c * lookupN.b;
    lk.Ap = R_c^2 * lookupN.Ap;
    lk.perim = R_c * lookupN.perim;
    lk.dAp_db = lk.perim;                       % dAp/db = perimeter
    lk.dperim_db = lookupN.dperim_db;
    % C12: never burn closer than web_min to the casing. The thinnest web sits
    % over the star tip, at R_c*(1 - x2) from the initial surface.
    b_stop = R_c*(1 - cfg.x2) - cfg.web_min;
    lk.b_max = min(R_c * lookupN.b_max, b_stop);

    if ~(lk.b_max > lk.b(1))
        r = struct("ok", false, "err", "C12 leaves no burnable web.");
        return
    end

    sim.grain_length = L;
    sim.throat_area = At;
    sim.mdot_ox = mdot_ox;
    sim.a_rf = cfg.a;
    sim.n_rf = cfg.n;
    sim.rho_f = cfg.rho_f;
    sim.Tc_fun = ctx.thermo.Tc_fun;
    sim.R_fun = ctx.thermo.R_fun;
    sim.k_fun = ctx.thermo.k_fun;
    sim.dRTdp_fun = ctx.thermo.dRTdp_fun;
    sim.dRTdOF_fun = ctx.thermo.dRTdOF_fun;
    sim.of_min = K.C5.lo;
    sim.of_max = K.C5.hi;
    sim.p_min = cfg.p_min;
    sim.p_max = cfg.p_max;
    sim.time_output = cfg.time_output;
    sim.fine_ode = cfg.fine_ode;
    sim.eps = ctx.engine.eps;
    sim.pamb = ctx.engine.pamb;

    r = run_combustion(lk, sim);
end
