# Risultati dell'ottimizzazione

Generato da `Optimization2/main_optimization.m` il 2026-08-20 18:32.
Non modificare a mano: viene riscritto a ogni run.

## 1. Vincoli attivi

Sorgente unica: `Optimization2/general/optimizationConstraints.m`. Nessun altro file contiene una soglia numerica di vincolo. Nessuna penalita' pesata: ogni vincolo e' un bound o un gate booleano.

| id | fase | tipo | grandezza | limite | motivazione | chi lo consuma, e contro cosa tira |
|----|------|------|-----------|--------|-------------|-------------------------------------|
| **C1** | A | bound | `h = x2 - x1/cos(pi/N)` | `>= 0` | Geometric: below h = 0 the tip radius falls inside the inner polygon vertex and the star does not exist. Reparametrizing (x1, x2) as (x1, h) turns the old linear inequality into this plain bound, so the optimum no longer sits on the face of a constraint. | Consumed by the phase A search box on h. Superseded in practice by C2, which floors h strictly above zero, so C1 never binds on its own. |
| **C2** | A | bound | `h / h_grid` | `>= 5` cells | The tip must be at least 5 grid cells tall, h_grid = R_c/grid_divisions. In the current design the optimal tip is 1/36 of a cell: what is being optimized there is discretization noise, not geometry. | h_min = C2.lo/grid_divisions, the lower bound of h in the phase A box. DEPENDS ON THE SEARCH GRID, NOT ON THE RE-RANKING GRID: with grid_search = 350 the floor is 0.01429, with grid_fine = 700 it is 0.00714, and the objective wants h -> 0, so the optimum sits exactly on whichever floor is in force. The shapes must therefore be RE-OPTIMIZED on the fine grid, not re-scored: a re-scored star would be ranked at twice the floor the ranking declares. This is also asymmetric between the shape families, because the cylinder has no h and no grid-derived bound at all: refining the grid moves the star and leaves the cylinder still, which biases the comparison in favour of the cylinder. C2 pulls against the objective, and against the fairness of the star-versus-cylinder comparison. |
| **C3** | A | bound | `max(x2, x1/cos(pi/N))` | `<= 1` | Geometric: the grain cannot stick out of the casing. | Consumed as the single linear inequality of the phase A search, x1/cos(pi/N) + h <= 1. Pulls against C2: a taller tip forces a smaller port, and both cannot grow together. |
| **C4** | A | gate | `G_ox(0) in phase A` | `[420, 666.667]` kg/(m2 s) | C10 applied already in phase A, on the flow the thrust requirement implies. Substituting R_c from C6 into G_ox(0) = mdot_ox/(R_c^2 Ap~(0)) factorizes it exactly as Gox0 = C(mdot_ox)*Gtilde, with C = mdot_ox^(1/(2n+1))/(t_b a)^(2/(2n+1)) carrying all the scale and Gtilde = I~(b~_end)^(2/(2n+1))/Ap~(0) depending only on the normalized shape; with n = 0.75 the exponents are 0.4 and 0.8. mdot_ox is NOT a free variable: C7 nails it down through mdot_ox = F_target/mean_t[(1 + 1/O_F)*c_eff(O_F)], every term of which phase A already has. DO NOT go back to the existential form, 'does some mdot_ox in a search box put G_ox(0) in band': it approves shapes that only work at flows the thrust forbids, and phase A parks on that corner. Observed symptom: x1 stuck at 0.3468 and G_ox(0) = 255 against a floor of 400. The band is C10's, tightened by the margin below, because phase A runs at a constant p_c = 20 bar and so gets mdot_ox wrong by a few percent. Derived, see below. | Consumed inside shapeMerit, level by level: an O/F level whose implied flow puts G_ox(0) out of band is discarded, and a shape survives only if some level does. THIS IS WHAT DECIDES THE PORT FRACTION x1: the objective is monotone in port size, so x1 has no interior optimum and ends up wherever C4 stops it. Pulls directly against the objective, and depends on C7 (which fixes the flow) and on C6 (which fixes R_c). |
| **C5** | A+B | gate | `O/F(t), whole history` | `[1.1, 19.5]` | Validity of the thermochemical lookup. Enforced over the WHOLE burn, not just at t = 0: outside this range the interpolants are extrapolating and the performance numbers mean nothing. | Consumed twice: as a truncation of the burn (the upper crossing is one of the candidate ends of burn, together with C12 and burnout) and as the clamp on every thermochemical call. Bounds the lambda sweep of phase A from both sides, since the initial O/F must land inside it. |
| **C6** | B | equation | `t_b` | `= 300` s | Assignment, not negotiable. Solved for, not searched: it determines R_c through t_b = R_c^(2n+1) I~(b~_end)/(a mdot_ox^n). | Consumed twice: as a closed form in phase A, where it gives R_c the moment mdot_ox is known, and as a bisection on R_c in sizeEngine. Feeds C4 (through the scale factor) and C12 (which needs R_c to turn 3 mm into a normalized length). Depends on b~_end, which depends on C5 and C12, so the three are resolved together. |
| **C7** | B | equation | `F` | `= 50000` N | Assignment. Determines the grain length L. Point (i) reads the 50 kN as the INITIAL thrust; the updated note allows reading them as the MEAN thrust when optimizing, in which case C9 becomes mandatory. Which reading is in force is the 'value' field of this row: 'mean' or 'initial'. It is the only switch in the table. | Consumed twice: in phase A it FIXES THE OXIDIZER FLOW through mdot_ox = F/mean[(1 + 1/OF) c_eff], which is what makes C4 evaluable there; in phase B it is solved for L on a bracket. Pulls against C9: on the mean reading the initial thrust drifts upwards and can leave the 45-55 kN band. On the mean reading it also fixes the total impulse at F*t_b = 15 MNs regardless of shape, so maximizing Isp on the loaded mass becomes exactly minimizing that mass. |
| **C8** | B | equation | `max p_c` | `= 2e+06` Pa | Assignment: it is a ceiling, and the design leans on it because Isp grows with p_c. Determines the throat area, A_t = max_b[mdot(b) cstar(b)]/p_target. | Consumed in sizeEngine, where p_c against A_t has an exact slope of -1 and the update is a single division. Phase A does not solve it: it assumes p_c = 20 bar throughout, which is the approximation the 5 % margin of C4 exists to cover. |
| **C9** | B | gate | `F(0)` | `[45000, 55000]` N | Assignment. Mandatory whenever the 50 kN of C7 are read as a mean thrust rather than an initial one. | Consumed by feas_check in phase B only: phase A never checks it, which is one of the ways a shape can pass A and fail B. Pulls against C7 on the mean reading, and tightens as the O/F drift grows, since a larger drift means a larger spread between F(0) and F mean. |
| **C10** | B | gate | `G_ox(0)` | `[400, 700]` kg/(m2 s) | THIS IS THE CONSTRAINT THAT CLOSES THE PROBLEM. Lower end: the value suggested by the assignment (500) and the validity of the regression correlation at start-up. Upper end: blow-off. Without the floor the optimum runs away to a 1 cm fuel annulus in a 3 m casing (G_ox(0) ~ 58). With it, G_ox(0) fixes r_0, R_c comes out of the burn time, and the geometry closes on its own. | Consumed by feas_check in phase B, and anticipated in phase A as C4. Its floor ALWAYS BINDS, so it is what selects the port fraction; that is why the choice of 400 is justified with a curve (sweepGoxFloor) rather than with a weight. Pulls against the objective, which would rather have a larger port, and against C7, which caps the flow. |
| **C11** | B | gate | `G_ox(t), whole history` | `>= 10` kg/(m2 s) | Not physics: rf = a G_ox^0.75 -> 0 stalls the ODE and the burn never terminates. A pure guard on the integrator. | Consumed by feas_check. Never binds once C10 holds, since the flux ends around 35 and the floor is 10. Kept as a guard, not as a design constraint. |
| **C12** | B | bound | `web_residual` | `>= 0.003` m | The grain is not burnt until the casing is exposed: liner plus insulation have to survive. Without it the objective rewards a zero sliver that cannot be built. | Consumed as an END OF BURN, in both phases: b~_end is the first of the C5 crossing, this stop, and burnout. Needs R_c from C6 to become a normalized length, so the two are iterated together. It is what makes the sliver comparable between shape families, and therefore what carries the cylinder-versus-star margin: the star's tips reach the casing before its flats, so at a common 3 mm of thinnest web the star still has fuel left elsewhere. |
| **C13** | B | gate | `2*R_c` | `<= 2` m | Safety net only. At G_ox(0) >= 200 the casing is 0.62 m across and about 6 m long, so this never binds, but it costs nothing to keep. ASSUMPTION: the assignment gives no envelope; these two numbers are ours, sized to leave a wide margin on the expected design (~0.60 m x ~6.1 m). Change them here and nowhere else. | Consumed by feas_check. Does not bind on HTPB with an O/F near 2, but it is the constraint most likely to bite on the high-O/F oxidizers, whose grains are much shorter but whose casings are wider. If it ever binds, the oxidizer ranking can invert, and that is a result to show rather than a failure. |
| | | | `L` | `<= 12` m | | |
| **C14** | - | removed | `oxidizer` | **RITIRATO** | REMOVED, and kept in the table so the reason is on record. 'GOX' in the assignment always denotes the oxidizer mass FLUX G_ox [kg/(m2 s)], not gaseous oxygen: it appears as 'partire con GOX = 500 kg/m2s', which is a flux. The oxidizer is therefore a free parameter at point (i) too, 'nel limite della ragionevolezza', and the best one becomes the baseline directly. Symmetry argument: the geometry is not imposed at point (i) either, so neither is the oxidizer; the levers of point (iii) are changes with respect to the baseline the student picks. GASEOUS oxygen is nonetheless excluded on storage grounds: ~3 t of oxidizer in the gas phase at 200 bar needs over 10 m3 of high-pressure tankage, which is a test-bench oxidizer, not an upper stage. Use O2(L). Justifying that exclusion is part of the critical analysis the assignment asks for. | Nothing consumes it any more. The oxidizer is now enumerated by main_optimization over the WHOLE A -> B chain, because it moves the peak of Isp(O/F), hence the best O/F level, hence mdot_ox, hence G_ox(0), hence the window C4 leaves open on the port fraction. |

**Derivati alle impostazioni correnti.** C2 dipende dalla griglia, e le fasi ne usano due: ricerca `grid = 350` -> `h >= 0.01429 R_c`, classifica `grid = 700` -> `h >= 0.00714 R_c`. Le forme sono **riottimizzate** sulla griglia di classifica, non solo rivalutate: il floor si dimezza fra le due e l'ottimo ci sta esattamente sopra. C4: `G_ox(0) in [420, 666.667]`, cioe' C10 con margine 5%, valutato sulla portata calcolata da C7. C7: i 50 kN letti come spinta **MEAN**.

## 2. Classifica delle forme, per ossidante

Fase A, a `R_c = 1`. `Isp_load` e' calcolata direttamente come `I_tot/(g0*m_load)`; `mdot_ox` non e' cercata, e' ricavata dal requisito di spinta C7. Il cilindro compete alla pari con le stelle ed e' valutato con la stessa lookup MDF.

### O2(L)

| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |
|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|
| 1 | cylinder | 0 | 0.2843 | 0.00000 | 337.22 | 2.043 | 0.0214 | 1.866 | inf | 10.08 | 422 |
| 2 | star | 18 | 0.2761 | 0.00714 | 336.59 | 2.047 | 0.0262 | 1.888 | inf | 10.09 | 423 |
| 3 | star | 17 | 0.2767 | 0.00714 | 336.53 | 2.055 | 0.0269 | 1.883 | inf | 10.10 | 420 |
| 4 | star | 16 | 0.2707 | 0.00714 | 336.46 | 2.048 | 0.0269 | 1.900 | inf | 10.09 | 437 |
| 5 | star | 15 | 0.2746 | 0.00714 | 336.42 | 2.048 | 0.0278 | 1.883 | inf | 10.09 | 422 |
| 6 | star | 14 | 0.2746 | 0.00714 | 336.38 | 2.053 | 0.0283 | 1.879 | inf | 10.10 | 420 |
| 7 | star | 13 | 0.2727 | 0.00714 | 336.28 | 2.049 | 0.0291 | 1.881 | inf | 10.09 | 423 |
| 8 | star | 12 | 0.2708 | 0.00714 | 336.20 | 2.049 | 0.0297 | 1.882 | inf | 10.09 | 425 |
| 9 | star | 11 | 0.2703 | 0.00714 | 336.15 | 2.050 | 0.0303 | 1.878 | inf | 10.09 | 421 |
| 10 | star | 10 | 0.2680 | 0.00714 | 336.04 | 2.050 | 0.0312 | 1.879 | inf | 10.09 | 422 |

### N2O

| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |
|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|
| 1 | cylinder | 0 | 0.3043 | 0.00000 | 304.14 | 5.728 | 0.0195 | 1.805 | inf | 14.23 | 421 |
| 2 | star | 18 | 0.2919 | 0.00714 | 303.84 | 5.739 | 0.0246 | 1.836 | inf | 14.24 | 434 |
| 3 | star | 17 | 0.2954 | 0.00714 | 303.83 | 5.739 | 0.0253 | 1.822 | inf | 14.24 | 421 |
| 4 | star | 16 | 0.2948 | 0.00714 | 303.81 | 5.739 | 0.0257 | 1.820 | inf | 14.24 | 421 |
| 5 | star | 15 | 0.2941 | 0.00714 | 303.78 | 5.742 | 0.0265 | 1.819 | inf | 14.24 | 421 |
| 6 | star | 14 | 0.2928 | 0.00714 | 303.76 | 5.743 | 0.0270 | 1.819 | inf | 14.24 | 422 |
| 7 | star | 13 | 0.2914 | 0.00714 | 303.73 | 5.743 | 0.0276 | 1.820 | inf | 14.24 | 423 |
| 8 | star | 12 | 0.2902 | 0.00714 | 303.70 | 5.744 | 0.0282 | 1.819 | inf | 14.24 | 423 |
| 9 | star | 11 | 0.2879 | 0.00714 | 303.66 | 5.745 | 0.0290 | 1.820 | inf | 14.24 | 425 |
| 10 | star | 10 | 0.2872 | 0.00714 | 303.62 | 5.747 | 0.0299 | 1.816 | inf | 14.24 | 420 |

### N2O4

| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |
|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|
| 1 | cylinder | 0 | 0.2946 | 0.00000 | 318.81 | 3.032 | 0.0204 | 1.834 | inf | 11.96 | 420 |
| 2 | star | 18 | 0.2858 | 0.00714 | 318.34 | 3.039 | 0.0254 | 1.855 | inf | 11.97 | 422 |
| 3 | star | 17 | 0.2852 | 0.00714 | 318.30 | 3.039 | 0.0259 | 1.854 | inf | 11.97 | 423 |
| 4 | star | 16 | 0.2816 | 0.00714 | 318.24 | 3.040 | 0.0262 | 1.863 | inf | 11.98 | 432 |
| 5 | star | 15 | 0.2832 | 0.00714 | 318.21 | 3.040 | 0.0270 | 1.854 | inf | 11.97 | 425 |
| 6 | star | 14 | 0.2836 | 0.00714 | 318.18 | 3.040 | 0.0276 | 1.849 | inf | 11.97 | 421 |
| 7 | star | 13 | 0.2825 | 0.00714 | 318.12 | 3.041 | 0.0283 | 1.848 | inf | 11.97 | 421 |
| 8 | star | 12 | 0.2802 | 0.00714 | 318.07 | 3.041 | 0.0288 | 1.851 | inf | 11.98 | 424 |
| 9 | star | 11 | 0.2799 | 0.00714 | 318.02 | 3.046 | 0.0296 | 1.846 | inf | 11.98 | 420 |
| 10 | star | 10 | 0.2766 | 0.00714 | 317.93 | 3.043 | 0.0305 | 1.850 | inf | 11.98 | 424 |

### IRFNA

| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |
|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|
| 1 | cylinder | 0 | 0.2997 | 0.00000 | 308.15 | 3.831 | 0.0201 | 1.819 | inf | 13.07 | 420 |
| 2 | star | 18 | 0.2906 | 0.00714 | 307.75 | 3.830 | 0.0251 | 1.840 | inf | 13.07 | 423 |
| 3 | star | 17 | 0.2909 | 0.00714 | 307.73 | 3.830 | 0.0256 | 1.836 | inf | 13.07 | 420 |
| 4 | star | 16 | 0.2904 | 0.00714 | 307.71 | 3.839 | 0.0260 | 1.834 | inf | 13.07 | 420 |
| 5 | star | 15 | 0.2897 | 0.00714 | 307.66 | 3.846 | 0.0268 | 1.833 | inf | 13.08 | 420 |
| 6 | star | 14 | 0.2887 | 0.00714 | 307.63 | 3.831 | 0.0274 | 1.832 | inf | 13.07 | 420 |
| 7 | star | 13 | 0.2874 | 0.00714 | 307.59 | 3.832 | 0.0279 | 1.832 | inf | 13.07 | 421 |
| 8 | star | 12 | 0.2862 | 0.00714 | 307.55 | 3.832 | 0.0284 | 1.831 | inf | 13.07 | 421 |
| 9 | star | 11 | 0.2841 | 0.00714 | 307.49 | 3.834 | 0.0293 | 1.833 | inf | 13.07 | 422 |
| 10 | star | 10 | 0.2813 | 0.00714 | 307.43 | 3.835 | 0.0301 | 1.835 | inf | 13.07 | 424 |

### H2O2_90

| # | forma | N | x1 | h | Isp_load [s] | O/F med | sigma | drift | drift a burnout | mdot_ox [kg/s] | G_ox(0) |
|---|-------|---|----|---|--------------|---------|-------|-------|-----------------|----------------|--------|
| 1 | cylinder | 0 | 0.2965 | 0.00000 | 344.90 | 5.405 | 0.0203 | 1.828 | inf | 12.43 | 421 |
| 2 | star | 18 | 0.2887 | 0.00714 | 344.58 | 5.417 | 0.0254 | 1.846 | inf | 12.44 | 420 |
| 3 | star | 17 | 0.2851 | 0.00714 | 344.55 | 5.418 | 0.0256 | 1.854 | inf | 12.44 | 429 |
| 4 | star | 16 | 0.2875 | 0.00714 | 344.53 | 5.420 | 0.0263 | 1.843 | inf | 12.44 | 420 |
| 5 | star | 15 | 0.2867 | 0.00714 | 344.50 | 5.421 | 0.0269 | 1.842 | inf | 12.44 | 421 |
| 6 | star | 14 | 0.2856 | 0.00714 | 344.47 | 5.422 | 0.0276 | 1.842 | inf | 12.44 | 421 |
| 7 | star | 13 | 0.2845 | 0.00714 | 344.44 | 5.424 | 0.0281 | 1.841 | inf | 12.44 | 421 |
| 8 | star | 12 | 0.2834 | 0.00714 | 344.40 | 5.426 | 0.0288 | 1.840 | inf | 12.44 | 421 |
| 9 | star | 11 | 0.2820 | 0.00714 | 344.35 | 5.428 | 0.0296 | 1.839 | inf | 12.44 | 420 |
| 10 | star | 10 | 0.2793 | 0.00714 | 344.31 | 5.429 | 0.0304 | 1.841 | inf | 12.44 | 422 |

## 3. Comparativa fra ossidanti

| ossidante | forma | Isp_load [s] | mdot_ox [kg/s] | O/F med | 2R_c [m] | L [m] | L/D | G_ox(0) | G_ox(fine) | sigma | esito |
|-----------|-------|--------------|----------------|---------|----------|-------|-----|---------|------------|-------|-------|
| O2(L) | cylinder | **336.38** | 10.08 | 2.029 | 0.613 | 6.097 | 9.95 | 422 | 34.8 | 0.0212 | OK |
| N2O | cylinder | **303.86** | 14.27 | 5.733 | 0.683 | 2.496 | 3.66 | 421 | 39.7 | 0.0193 | OK |
| N2O4 | cylinder | **318.29** | 11.96 | 3.006 | 0.646 | 4.422 | 6.84 | 420 | 37.1 | 0.0202 | OK |
| IRFNA | cylinder | **307.79** | 13.07 | 3.829 | 0.664 | 3.602 | 5.42 | 420 | 38.4 | 0.0198 | OK |
| H2O2_90 | cylinder | **344.53** | 12.43 | 5.302 | 0.654 | 2.548 | 3.89 | 421 | 37.7 | 0.0200 | OK |

Ossidanti enumerati: O2(L), N2O, N2O4, IRFNA, H2O2_90.
Un ossidante che non dimensiona e' riportato con il motivo, mai scartato in silenzio.

## 4. Perche' ha vinto

Confronto in Fase A per **H2O2_90**, migliore cilindro contro migliore stella, valutati con la stessa lookup MDF e la stessa griglia.

| | cilindro | stella (N = 18) |
|---|---------|--------|
| `x1` | 0.2965 | 0.2887 |
| `Isp_load` [s] | **344.90** | 344.58 |
| `Isp_med` [s] (drift) | 346.01 | 345.98 |
| fattore sliver | 0.99678 | 0.99595 |
| `sigma` | 0.0203 | 0.0254 |
| drift (parte utile) | 1.828 | 1.846 |
| drift fino a burnout | inf | inf |

**Scomposizione della differenza di +0.32 s:**

- da drift (via `Isp_med`): **+0.03 s**
- da sliver (via il fattore di massa): **+0.29 s**
- somma: +0.32 s

Il termine dominante e' lo **sliver**.

**Sul drift fino a burnout.** La colonna diverge per entrambe le famiglie, cilindro compreso: il solutore MDF smette di contare la superficie che esce dal casing, quindi il perimetro bruciante va a zero alla parete per qualunque forma e `Phi = Ap^n/Pb` diverge con lui. Non e' un fenomeno fisico della stella, e' contabilita' di fine burn: le zone di collasso sono confrontabili fra le due forme (ordine dell'1% del web). Il burn reale si ferma prima, su C12, e li' il drift e' quello della colonna 'parte utile'. La colonna a burnout va riportata al punto (ii) come limite del modello, non come proprieta' della geometria.

**Quello che invece distingue davvero le due forme e' lo sliver al criterio di arresto C12.** Le punte della stella raggiungono il casing prima dei piani, quindi quando il punto piu' sottile del web e' sceso a 3 mm resta ancora combustibile altrove. Il cilindro e' l'unica porta concentrica che raggiunge la parete ovunque nello stesso istante, e questo si vede direttamente nel confronto di `sigma`: 0.0203 contro 0.0254.

---

*Decomposizione esatta:* `Isp_load = Isp_med * (OF_med + 1)/(OF_med + 1/(1 - sigma))`. Il primo fattore misura la perdita da drift dell'O/F, il secondo quella da sliver. E' una diagnostica: `Isp_load` viene calcolata direttamente da `I_tot/(g0*m_load)`, non per questa strada.
