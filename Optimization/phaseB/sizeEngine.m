function S = sizeEngine(mdot_ox, ctx, lookupN, cfg)
    % sizeEngine
    % Phase B sizing: given the oxidizer flow and the normalized grain shape,
    % solve the three equations C6, C7, C8 for the three sizes they determine.
    %
    %   C6  burn time     t_b = 300 s      ->  casing radius R_c
    %   C7  thrust        F   = 50 kN      ->  grain length L
    %   C8  peak pressure p_c = 20 bar     ->  throat area A_t
    %
    % Each one is solved on a bracket with a verified sign change, never by a
    % damped iteration: the end of burn b~_end itself depends on R_c and on L
    % through the C5 and C12 events, so t_b(R_c) is not guaranteed monotone and
    % a damped update can walk off. The autonomous-burnback quadrature
    %   t(b~) = R_c^(2n+1) I~(b~)/(a mdot_ox^n)
    % is used only for the initial guesses and the brackets; every value that
    % gets reported comes from run_combustion, the same exact ODE model the
    % final validation run uses.
    %
    % C12 is applied as a hard stop on the burnback: the grain is never burnt
    % closer than web_min to the casing, so the event fires at
    %   b_stop = R_c*(1 - x2) - web_min
    % which is the fuel thickness over the casing at the star tip, the thinnest
    % point of the web.
    % INPUT
    %   mdot_ox / double / 1x1   Constant oxidizer mass flow [kg/s]
    %   ctx     / struct / 1x1   Context from engineContext (K, thermo, engine)
    %   lookupN / struct / 1x1   Normalized geometry lookup at R_c = 1
    %   cfg     / struct / 1x1   n, a, rho_f [-, m/s, kg/m3], x2 [-] tip radius
    %                            ratio, web_min [m], p_min, p_max [Pa],
    %                            time_output [s], fine_ode, and optionally
    %                            max_iter, tol
    % OUTPUT
    %   S       / struct / 1x1   ok, err, the sized engine (R_c, L, At,
    %                            throat_diameter), its performance and the
    %                            diagnostics listed in the brief

    S = struct("ok", false, "err", "");

    K = ctx.K;
    n = cfg.n;
    a = cfg.a;
    t_target = K.C6.lo;                 % [s]   C6
    F_target = K.C7.lo;                 % [N]   C7
    p_target = K.C8.lo;                 % [Pa]  C8
    thrust_mode = K.C7.value;           % "mean" or "initial"

    if ~isfield(cfg, "max_iter") || isempty(cfg.max_iter), cfg.max_iter = 4; end
    if ~isfield(cfg, "tol") || isempty(cfg.tol), cfg.tol = 3e-3; end

    % ------------ INITIAL GUESSES FROM THE QUADRATURE ------------
    % R_c from C6 with b~_end at the full web, L from the mean shape function
    % at the best-Isp mixture ratio, A_t from the choked-throat relation
    I_tilde = trapz(lookupN.b, lookupN.Ap.^n);              % [-]  burn integral
    if ~(I_tilde > 0)
        S.err = "Degenerate burn integral on the normalized lookup.";
        return
    end
    R_c = (t_target * a * mdot_ox^n / I_tilde)^(1/(2*n + 1));       % [m]

    % Mean shape function as a RATIO OF INTEGRALS, not as the mean of the
    % ratio: Phi~ = Ap~^n/Pb~ diverges at the end of the burn, where the
    % perimeter collapses, so its arithmetic mean is meaningless. This form is
    % the one that reproduces the mass-weighted O/F exactly,
    %   OF_med = lambda * int(Ap~^n db~)/int(Pb~ db~)
    Phi_mean = I_tilde / max(trapz(lookupN.b, lookupN.perim), realmin);
    L = mdot_ox^(1 - n) * R_c^(2*n - 1) * Phi_mean / (cfg.rho_f * a * ctx.of_best);
    At = mdot_ox * (1 + 1/ctx.of_best) * ctx.cstar_best / p_target;   % [m2]

    % ------------ FIXED POINT ON THE THREE EQUATIONS ------------
    r = struct("ok", false);
    last_fail = "";                 % the model's own reason for the last refusal
    converged = false;
    for iter = 1:cfg.max_iter
        % C6 -> R_c, by bisection on a verified bracket
        [R_c, okR, msg] = solve_bracketed(@(x) burn_time_of(x, L, At), ...
            R_c, t_target, cfg.tol);
        if ~okR
            S.err = "C6 (burn time -> R_c): " + msg + " Model says: " + last_fail;
            return
        end

        % C7 -> L, on a verified bracket
        [L, okL, msg] = solve_bracketed(@(x) thrust_of(R_c, x, At), ...
            L, F_target, cfg.tol);
        if ~okL
            S.err = "C7 (thrust -> L): " + msg + " Model says: " + last_fail;
            return
        end

        % C8 -> At. p_c is inversely proportional to A_t at fixed mass flow, so
        % the ratio update is an exact Newton step: it converges in one or two
        % passes and needs no bracket.
        for k = 1:6
            r = simulate(R_c, L, At);
            if ~r.ok
                S.err = "C8 (peak pressure -> At): " + reason(r);
                return
            end
            ratio = r.max_pressure / p_target;
            if abs(ratio - 1) < cfg.tol
                break
            end
            At = At * ratio;
        end

        % Residuals of all three equations at the current point
        r = simulate(R_c, L, At);
        if ~r.ok
            S.err = "Final evaluation: " + reason(r);
            return
        end
        err_all = max([abs(r.burn_time/t_target - 1), ...
                       abs(thrust_value(r)/F_target - 1), ...
                       abs(r.max_pressure/p_target - 1)]);
        if err_all < cfg.tol
            converged = true;
            break
        end
    end

    if ~converged
        S.err = sprintf("Fixed point did not converge in %d passes (residual %.3g).", ...
            cfg.max_iter, err_all);
        return
    end

    % ------------ THE SIZED ENGINE ------------
    S.ok = true;
    S.mdot_ox = mdot_ox;
    S.R_c = R_c;                                  % [m]
    S.L = L;                                      % [m]
    S.At = At;                                    % [m2]
    S.throat_diameter = sqrt(4*At/pi);            % [m]
    S.b_stop = R_c*(1 - cfg.x2) - cfg.web_min;    % [m]  C12 limit

    S.burn_time = r.burn_time;                    % [s]
    S.mean_thrust = r.mean_thrust;                % [N]
    S.thrust0 = r.initial_thrust;                 % [N]
    S.p_peak = r.max_pressure;                    % [Pa]
    S.I_tot = r.total_impulse;                    % [Ns]
    S.m_prop = r.m_prop;                          % [kg]  through the nozzle
    S.of0 = r.of0; S.of_end = r.of_end;           % [-]
    S.Gox0 = r.Gox0; S.Gox_end = r.Gox_end;       % [kg/m2 s]

    % ------------ OBJECTIVE: Isp ON THE LOADED MASS ------------
    % The denominator is the mass that leaves the ground, not the mass that
    % crosses the nozzle: the fuel left on board at burnout is dead weight for
    % an upper stage and has to be paid for.
    g0 = 9.80665;
    Ap0 = R_c^2 * lookupN.Ap(1);                              % [m2]
    b_end = r.in_time.burnback(end);                          % [m]
    Ap_end = interp1(R_c*lookupN.b, R_c^2*lookupN.Ap, ...
        min(b_end, R_c*lookupN.b(end)), "linear");            % [m2]
    A_case = pi * R_c^2;                                      % [m2]

    m_ox_load = mdot_ox * r.burn_time;                        % [kg]
    m_f_load = cfg.rho_f * L * (A_case - Ap0);                % [kg]
    S.m_load = m_ox_load + m_f_load;                          % [kg]
    S.Isp_load = r.total_impulse / (g0 * S.m_load);           % [s]
    S.Isp_med = r.mean_Isp;                                   % [s]  burnt mass
    S.sigma = (A_case - Ap_end) / max(A_case - Ap0, realmin); % [-]  sliver
    S.OF_med = m_ox_load / max(cfg.rho_f*L*(Ap_end - Ap0), realmin);

    % ------------ DIAGNOSTICS (NOT constraints) ------------
    S.Gox_min = min(r.in_time.GOX);                           % [kg/m2 s]
    S.Gox_150 = interp1(r.in_time.time, r.in_time.GOX, 150, "linear", "extrap");
    S.drift = max(r.in_time.O_F) / max(min(r.in_time.O_F), realmin);
    S.of_min_hist = min(r.in_time.O_F);
    S.of_max_hist = max(r.in_time.O_F);
    S.rf0 = cfg.a * r.Gox0^n;                                 % [m/s]
    S.rf_end = cfg.a * r.Gox_end^n;                           % [m/s]
    S.web_residual = R_c*(1 - cfg.x2) - b_end;                % [m]  C12
    S.r = r;

    %% Nested helpers, closing over the fixed data of this sizing

    function t = burn_time_of(R_c_try, L_try, At_try)
        rr = simulate(R_c_try, L_try, At_try);
        if rr.ok
            t = rr.burn_time;
        else
            last_fail = reason(rr);
            t = NaN;
        end
    end

    function F = thrust_of(R_c_try, L_try, At_try)
        rr = simulate(R_c_try, L_try, At_try);
        if rr.ok
            F = thrust_value(rr);
        else
            last_fail = reason(rr);
            F = NaN;
        end
    end

    function F = thrust_value(rr)
        % C7 is read either as the mean thrust or as the initial one; which is
        % in force is the 'value' field of the C7 row, nowhere else
        if thrust_mode == "initial"
            F = rr.initial_thrust;
        else
            F = rr.mean_thrust;
        end
    end

    function rr = simulate(R_c_try, L_try, At_try)
        % The exact ODE model, through the shared evaluator so that the solver
        % and verifySolver cannot drift apart
        rr = simulateEngine(R_c_try, L_try, At_try, mdot_ox, ctx, lookupN, cfg);
    end
end

%% FUNCTIONS

function msg = reason(r)
    % reason
    % Extract the model's own explanation of a failed run, so a solver failure
    % is never reported as an infeasible design.
    % INPUT
    %   r   / struct / 1x1   run_combustion output
    % OUTPUT
    %   msg / string / 1x1   Explanation, or a generic message

    if isfield(r, "err") && strlength(string(r.err)) > 0
        msg = string(r.err);
    else
        msg = "run_combustion returned no result and no reason.";
    end
end

function [x, ok, msg] = solve_bracketed(fun, x0, target, tol_rel)
    % solve_bracketed
    % Solve fun(x) = target on a bracket whose sign change has been verified,
    % expanding the bracket geometrically around the initial guess until one is
    % found. Bracketed by construction: the iterate can never leave the
    % interval, which a damped fixed point can.
    % INPUT
    %   fun     / handle / 1x1   Monotone-ish target function of one variable
    %   x0      / double / 1x1   Initial guess (must be positive)
    %   target  / double / 1x1   Value to hit
    %   tol_rel / double / 1x1   Relative tolerance on the residual [-]
    % OUTPUT
    %   x       / double / 1x1   Solution
    %   ok      / logical/ 1x1   True when a bracket was found and solved
    %   msg     / string / 1x1   Why it failed, empty on success

    ok = false;
    msg = "";
    x = x0;
    g = @(v) fun(v) - target;

    % The guess may land where the model refuses to run at all. Probe around it
    % before giving up: a bad guess is not the same thing as a bad design.
    g0 = g(x0);
    if ~isfinite(g0)
        probes = [1.02, 1/1.02, 1.05, 1/1.05, 1.15, 1/1.15, 1.4, 1/1.4, ...
                  2, 1/2, 4, 1/4, 10, 1/10];
        for p = probes
            x_try = x0 * p;
            g_try = g(x_try);
            if isfinite(g_try)
                x0 = x_try;
                g0 = g_try;
                break
            end
        end
    end
    if ~isfinite(g0)
        msg = sprintf(['the model does not evaluate anywhere within a factor 10 ' ...
            'of the initial guess %.4g.'], x);
        return
    end
    x = x0;
    if abs(g0) <= tol_rel*abs(target)
        ok = true;
        return
    end

    % Walk outwards from the guess until the residual changes sign. The step
    % ACCELERATES while the model keeps evaluating and BACKS OFF when it stops,
    % instead of following a fixed ladder of factors. The difference matters:
    % the model has a hard feasibility edge (the O/F leaves the CEA table a few
    % percent away from the design point), and a first probe at +30 % lands
    % beyond it, reports NaN, and hides a root sitting at +1 %.
    going_up = g0 < 0;                   % undershooting: the root is above
    x_good = x0;
    g_good = g0;
    step = 1.02;                         % [-] initial relative step
    found = false;
    for k = 1:60
        if going_up
            x_try = x_good * step;
        else
            x_try = x_good / step;
        end
        g_try = g(x_try);

        if isfinite(g_try)
            if sign(g_try) ~= sign(g0)
                % Bracket closed between the last good point and this one
                if going_up
                    lo = x_good; glo = g_good; hi = x_try; ghi = g_try;
                else
                    lo = x_try; glo = g_try; hi = x_good; ghi = g_good;
                end
                found = true;
                break
            end
            x_good = x_try;
            g_good = g_try;
            step = 1 + min((step - 1)*1.6, 0.5);      % accelerate, cap at +50 %
        else
            step = 1 + (step - 1)/3;                  % back off towards the edge
            if step < 1.0005
                break                                  % the edge is here
            end
        end

        if x_good > 1e6*x0 || x_good < 1e-6*x0
            break
        end
    end

    if ~found
        msg = sprintf(['no sign change between %.4g and %.4g (residual %.4g at ' ...
            'the far end): the target is out of reach on this side, the model ' ...
            'stops evaluating first.'], x0, x_good, g_good);
        return
    end

    % Bisection on the verified bracket
    for k = 1:40
        mid = 0.5*(lo + hi);
        gm = g(mid);
        if ~isfinite(gm)
            msg = sprintf("the model stopped evaluating inside the bracket, at %.4g.", mid);
            return
        end
        if abs(gm) <= tol_rel*abs(target) || (hi - lo) <= 1e-6*mid
            x = mid;
            ok = true;
            return
        end
        if gm > 0
            hi = mid;
        else
            lo = mid;
        end
    end

    x = 0.5*(lo + hi);
    ok = true;      % bracket never left, residual simply not tightened further
end
