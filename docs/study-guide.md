# Understanding and defending the SDOF project

The objective is to explain the modeling decisions, derive the response, modify the implementation, and judge the evidence. Memorizing a damping ratio or reading code aloud is insufficient. Work through the exercises before reading their checks. Use a separate scratch script for experiments so the published case study stays reproducible.

## 1. Physical model: what system are you actually studying?

Imagine a mass on a frictionless guide, attached to a fixed support by a spring and a viscous damper in parallel. Pull it away from equilibrium, hold it at rest, and release it. The spring accelerates the mass toward equilibrium; inertia carries it onward; the damper dissipates energy whenever the mass moves.

The model has one degree of freedom because one coordinate, $x(t)$, specifies its configuration. Velocity is a state variable, not another geometric degree of freedom. An SDOF second-order system therefore has two states.

Choose positive displacement to the right. The spring force is $-kx$. The damping force is $-c\dot x$. The damper opposes velocity, which need not have the same sign as displacement. Newton's law gives

$$m\ddot x=-c\dot x-kx.$$

Every term has units of force:

$$[m\ddot x]=\mathrm{kg\,m/s^2},\quad
[c\dot x]=(\mathrm{N\,s/m})(\mathrm{m/s}),\quad[kx]=(\mathrm{N/m})(\mathrm m).$$

If vertical displacement $q$ is measured downward from the unloaded spring length, the equation includes $mg$. Set $q=q_{\rm eq}+x$ and $kq_{\rm eq}=mg$; the constant gravity term cancels. Gravity has not disappeared physically; the coordinate origin incorporates static equilibrium.

**Be able to state the assumptions:** fixed support; one translational coordinate; constant positive mass and stiffness; constant nonnegative viscous damping; linear restoring force; a displacement range over which linearity is accepted; no external time-varying input after release. The given 0.50 m displacement is illustrative. Without an actual geometry you cannot demonstrate that it is small enough for a physical device.

**Exercise:** draw the free-body diagram for $x>0,\dot x<0$. Which way does each force point? Then derive the equation for a vertical system using $q$ before introducing $x$.

**Check:** the spring points left; the damper points right. The net force depends on their magnitudes.

## 2. Parameters and nondimensionalization

Divide by mass:

$$\ddot x+\frac{c}{m}\dot x+\frac{k}{m}x=0.$$

Define

$$\omega_n=\sqrt{k/m},\quad c_{\rm crit}=2\sqrt{km},\quad
\zeta=\frac{c}{c_{\rm crit}},\quad \frac{c}{m}=2\zeta\omega_n.$$

$\omega_n$ is an angular frequency in rad/s; $f_n=\omega_n/(2\pi)$ is cycles/s. The coefficient $c$ is dimensional. The damping ratio $\zeta$ compares the actual damping to the critical value for the chosen mass and stiffness. It is not the fraction of energy lost per cycle.

For $x_0\ne0$, let $y=x/x_0$ and $\tau=\omega_n t$. The chain rule gives

$$\dot x=x_0\omega_n y',\quad\ddot x=x_0\omega_n^2y'',\quad
y''+2\zeta y'+y=0.$$

The normalized initial conditions are

$$y(0)=1,\quad y'(0)=\frac{v_0}{\omega_n x_0}.$$

Thus response shape depends on $\zeta$ **and normalized initial velocity**. For release from rest, normalized initial velocity is fixed at zero. Frequency sets the dimensional time scale.

**Prediction exercises:**

1. Double $x_0$ with $v_0=0$. Predict displacement, normalized displacement, relative-band settling time, and energy.
2. Multiply $m,k,c$ by the same positive factor. Does the displacement history change?
3. Double $k$ while keeping $m,c$ fixed. What happens to $\omega_n$ and $\zeta$?

**Checks:** (1) displacement doubles, normalized response and relative settling time stay fixed, energy quadruples. (2) The ODE after division by mass is unchanged. (3) $\omega_n$ increases by $\sqrt2$ and $\zeta$ decreases by $\sqrt2$; this is not an experiment at constant damping ratio.

## 3. Derive the three response regimes

Try $x=e^{rt}$ because derivatives of an exponential are multiples of the same exponential. Substitution leads to

$$r^2+2\zeta\omega_nr+\omega_n^2=0,$$

$$r_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.$$

### Underdamped: $0\leq\zeta<1$

Set $\alpha=\zeta\omega_n$ and $\omega_d=\omega_n\sqrt{1-\zeta^2}$. Euler's identity converts conjugate complex exponentials into a real solution:

$$x=e^{-\alpha t}[A\cos(\omega_dt)+B\sin(\omega_dt)].$$

Apply $x(0)=x_0$ and differentiate once:

$$A=x_0,\quad v_0=-\alpha A+\omega_d B,\quad
B=\frac{v_0+\alpha x_0}{\omega_d}.$$

A common mistake is using $x_0e^{-\alpha t}\cos(\omega_dt)$ for release from rest. Its initial velocity is $-\alpha x_0$, not zero. The sine term is essential.

### Critical: $\zeta=1$

The repeated root $r=-\omega_n$ yields two independent solutions $e^{-\omega_nt}$ and $te^{-\omega_nt}$:

$$x=[x_0+(v_0+\omega_nx_0)t]e^{-\omega_nt}.$$

Explain why two arbitrary constants are necessary: the second-order equation must accommodate two independent initial conditions. Repeating the same exponential twice would supply only one independent solution.

### Overdamped: $\zeta>1$

$$x=C_1e^{r_1t}+C_2e^{r_2t},\quad
C_1+C_2=x_0,\quad r_1C_1+r_2C_2=v_0.$$

Solve the two equations to recover the coefficients in the derivation. Both roots are negative, but the slower root dominates late motion. For very large damping,

$$r_{\rm slow}\approx-\frac{k}{c},\qquad
r_{\rm fast}\approx-\frac{c}{m}.$$

Excess damping can therefore slow recovery. A nonoscillatory system may still cross equilibrium once for some initial velocities. The monotonic statements in this study are tied to release from rest.

**Mastery task:** derive all three solutions on paper, apply the initial conditions, and explain what changes at $\zeta=1$. Repeat with $x_0=0$ and $v_0\ne0$.

## 4. Energy and initial acceleration

The mass stores kinetic energy and the spring stores potential energy:

$$E=\frac12m\dot x^2+\frac12kx^2.$$

Differentiate and use the ODE:

$$\dot E=\dot x(m\ddot x+kx)=-c\dot x^2.$$

This has stronger content than a nice-looking curve: negative damping would add energy, zero damping conserves energy, and positive damping never increases it. With $c>0$, instantaneous dissipation is zero whenever velocity is zero. The energy trace need not be a straight exponential on a logarithmic plot because kinetic and potential energy exchange during motion.

The initial conditions imply

$$\ddot x(0)=-\frac{cv_0+kx_0}{m}.$$

For this experiment $v_0=0$, so every case starts at $-20$ m/s² regardless of damping. Differentiate once more to find $x^{(3)}(0)$ and obtain

$$x(t)=x_0\left[1-\frac{\omega_n^2t^2}{2}
+\frac{\zeta\omega_n^3t^3}{3}+O(t^4)\right].$$

**Exercise:** derive the cubic coefficient without expanding the closed-form solution. Explain why damping first appears at cubic order in displacement for release from rest.

## 5. Metrics: distinguish the questions you are answering

| Metric | Question | Main qualification |
| --- | --- | --- |
| Zero-crossing time | When does the mass first pass equilibrium? | It may move far away again. |
| First reverse peak | How large is the first excursion on the other side? | The simple formula here requires $v_0=0$. |
| 2% settling time | When does the response enter and remain within the band? | A finite sampled record only estimates this. |
| Envelope half-life | When does the exponential factor halve? | It is not the first time displacement halves. |
| Energy | How much mechanical energy remains? | It is not displacement or a comfort criterion. |

For underdamping with $v_0=0$, differentiate the solution and simplify:

$$\dot x=-\frac{x_0\omega_n^2}{\omega_d}e^{-\zeta\omega_nt}\sin(\omega_dt).$$

Extrema occur at integer multiples of $\pi/\omega_d$. At the first reverse peak,

$$t_p=\pi/\omega_d,\quad x(t_p)/x_0=-e^{-\pi\zeta/\sqrt{1-\zeta^2}}.$$

Define $M_p=|x(t_p)|/|x_0|$. To solve $M_p=\delta$, take logarithms, set $L=-\ln\delta>0$, and rearrange:

$$L=\frac{\pi\zeta}{\sqrt{1-\zeta^2}},\quad
L^2(1-\zeta^2)=\pi^2\zeta^2,\quad
\zeta=\frac{L}{\sqrt{\pi^2+L^2}}.$$

For $\delta=0.02$, this gives $\zeta=0.779703$, $c=9862.553$ N·s/m, and the sampled nominal settling estimate is about 0.569603 s. The negative peak occurs later, around 0.7933 s, but touches rather than exceeds the permitted boundary. Equality is allowed.

**Why the sweep has jumps:** a smoothly shrinking late peak eventually stops violating the band. The last disqualifying excursion then changes to an earlier one. The minimum feasible time under a threshold definition can jump even while every displacement value changes continuously.

**Why the boundary needs caution:** slightly reduce $\zeta$ and the negative peak exceeds 2%; settling must wait until after that excursion. Parameter uncertainty, a nonzero release velocity, or a different band can materially change the result. Do not call this a universally optimal suspension or aircraft damping ratio.

The common estimate $t_s\approx4/(\zeta\omega_n)$ is an envelope-based rule, not the exact last crossing for this initial-value problem. For $v_0=0$, a conservative underdamped displacement envelope gives

$$t\geq\frac{\ln[1/(\delta\sqrt{1-\zeta^2})]}{\zeta\omega_n}.$$

This is sufficient for all later displacement to lie in the band, but generally conservative.

**Experiment:** compare $\zeta=0.77$, the exact design, $0.79$, $1$, and $1.3$. Predict reverse peak and settling before plotting. Repeat with a 5% band without silently changing the 2% design label.

## 6. MATLAB: understand the language used here

Read `free_response.m` first, then the settling helper, then `run_analysis.m`, and finally the tests. Functions implement reusable calculations; the analysis script chooses the experiment.

| Code | Meaning |
| --- | --- |
| `function [x,v,a,meta] = free_response(...)` | Named function with four outputs; fewer can be requested. |
| `[x,v,~,meta]` | Request outputs but discard acceleration. |
| `@(~,state) [...]` | Anonymous function; time argument exists but is unused in this autonomous ODE. |
| `(0:1e-4:4).'` | Samples from 0 to 4 s; nonconjugate transpose makes a column. |
| `xHistory(:,i)` | All time samples for case `i`. |
| `.*`, `./`, `.^` | Elementwise operations; `*`, `/`, `^` have matrix meanings. |
| `zeros(numel(t),numberOfCases)` | Preallocate a matrix: rows are time, columns are cases. |
| `meta.zeta` | Read a named field from the output structure. |
| `...` | Continue a statement on the next source line. |
| `NaN` | Missing/not-applicable numerical value; understand its context. |
| `assert(condition,message)` | Stop execution if a required check fails. |
| `mfilename`, `fileparts`, `fullfile` | Locate the script and build portable file paths. |
| `addpath` | Allow MATLAB to locate functions in the project directory. |

The `~` used for ignored outputs and the `~` logical NOT operator have different roles. A semicolon suppresses command-window output; it does not alter the numerical result.

**Exercise:** create `t` with five samples and manually predict the size of each intermediate array. Write a small loop that fills a displacement column. Then replace it with one vectorized call and compare.

**Common mistakes:** typing Markdown backslashes into function names; using `^2` on a displacement vector; reading a row as a case when rows are times; running a relative path from inside `tests` and accidentally requesting `tests/tests`; opening files directly inside an unextracted ZIP.

## 7. Numerical formulation: what does ode45 actually do?

Introduce $q=[x\;v]^T$:

$$\dot q=Aq,\quad
A=\begin{bmatrix}0&1\\-k/m&-c/m\end{bmatrix},\quad
q(0)=\begin{bmatrix}x_0\\v_0\end{bmatrix}.$$

The first row says displacement changes at the velocity. The second row is Newton's law divided by mass. `ode45` integrates these coupled first-order equations using an adaptive Runge–Kutta method. It does not symbolically solve the characteristic polynomial.

The requested output times in `tVerify` do not force the internal solver step size. Relative and absolute tolerances control estimated local numerical error; they do not certify physical accuracy or guarantee the final error equals the tolerance. Absolute tolerance matters near zero, where relative error is ill-conditioned.

The comparison takes absolute displacement differences at matching output times and reports their maximum. A difference around $10^{-11}$ m means the two computations agree closely for the chosen equation. It does not mean a physical displacement is known to picometers.

For constant $A$, the state also satisfies $q(t)=e^{At}q(0)$. MATLAB `expm(A*t)` provides another computational route, so tests compare both displacement and velocity for multiple initial states and damping ratios. It is still a solution of the same model, not an experiment.

**Exercise:** write the anonymous derivative function from memory. Explain why `[x0;v0]` is a column and why the second output column from `ode45` is velocity. Increase solver tolerances and observe how the discrepancy changes.

## 8. Floating-point robustness: why the code looks different from the textbook

Near critical damping, $\omega_d$ becomes small. Computing a huge coefficient $B$ and multiplying it by a tiny sine can lose numerical quality. Instead the implementation uses

$$\frac{\sin(\omega_dt)}{\omega_d}=t\frac{\sin u}{u},\quad u=\omega_dt,$$

and a series $\sin u/u\approx1-u^2/6+u^4/120$ near zero. This is a local numerical evaluation technique, not a truncation of the full physical solution over the entire time span. MATLAB's built-in `sinc` has a different $\pi$ normalization.

The overdamped implementation rationalizes the slow root to avoid subtracting two nearly equal positive magnitudes. It uses `expm1` to evaluate differences of nearby exponentials accurately. Derive the stable form in `derivation.md` from the usual two-exponential solution; do not memorize it blindly.

An eight-ULP neighborhood of $\zeta=1$ is treated as the repeated-root limit. ULP means the spacing between neighboring representable floating-point numbers. This tiny implementation convention should not be confused with the much larger engineering category “near critical.”

The settling helper also allows a tiny floating-point slack so evaluating the exact tangent peak does not falsely report a physically nonexistent excursion. The slack is about arithmetic, not permission to violate a design requirement by a meaningful amount.

**Exercise:** compare naive and stable root formulas at very large $\zeta$. Explain why the software explicitly does not promise reliable results for arbitrary inputs so extreme that derived quantities overflow.

## 9. Settling-time algorithm and its limitations

Read the helper and narrate it in order:

1. Validate matching nonempty vectors, at least two samples, increasing time, and a positive relative band.
2. Form the absolute band from the reference displacement.
3. Find the last sample outside the band, allowing only roundoff slack.
4. If none is outside, return the first supplied time; if the final sample is outside, return `NaN`.
5. Interpolate the crossing between the last outside sample and the next sample using the signed boundary.

For points $(t_1,x_1)$ and $(t_2,x_2)$ with boundary $b$,

$$t_s\approx t_1+\frac{b-x_1}{x_2-x_1}(t_2-t_1).$$

The function cannot see between samples or beyond the final time. A zero-valued sample is not proof of settling. Even all-inside samples could miss oscillations if sampled too coarsely. For the specific model, analytical decay bounds plus grid/horizon checks support the reported values.

**Exercise:** hand-evaluate the helper for `t=[0,1,2,3]`, `x=[1,0,0.1,0]`, 2% tolerance and reference 1. Then shorten the record to its first two samples. Explain why the second answer is misleading about the later full record.

**Check:** the full record gives 2.8 s; the short record estimates 0.98 s because it has no knowledge of re-excursion.

## 10. Verification: explain what each test can and cannot establish

| Check group | What failure would suggest |
| --- | --- |
| Initial conditions and initial acceleration | Incorrect coefficients, units, or state mapping; acceleration identity alone is not independent. |
| Analytical displacement vs ode45 | A branch, sign, parameter, or numerical integration mismatch. |
| Non-increasing energy | Wrong velocity, wrong damping sign, or numerical inconsistency. |
| Critical-limit continuity | Inconsistent branch limits or cancellation. |
| Nominal design regression | Changed convention or broken design/settling calculation. |
| Selected overdamped monotonicity | Sign reversal inconsistent with the selected release-from-rest case. |
| General initial states vs expm | Hidden reliance on $v_0=0$ or errors in either state. |
| Differential consistency and conservation | Derivative inconsistency or unphysical undamped energy change. |
| Settling edge cases | Incorrect last-crossing or invalid-input behavior. |
| Peak, grid, horizon, and tail checks | A sampling artifact, roundoff at tangency, or inadequate record duration. |

Tests are evidence within tested cases, not proof of all possible inputs. Checking $a=-(cv+kx)/m$ when that is exactly how `a` was computed is circular; this is why finite-difference derivative checks and independent state solutions matter.

**Verification:** did we implement and solve the chosen equations correctly?
**Validation:** do those equations adequately describe the physical system for the intended use?
This project contains verification. Validation awaits a physical device, measurements, calibrated parameters, and quantified uncertainty.

**Exercise:** intentionally flip the damper sign in a scratch copy. Predict which tests fail before running them. Restore the code afterward. Then halve the time-grid spacing and explain why unchanged plots alone are weaker evidence than numerical metric comparisons.

## 11. Repository and figure pipeline

`run_analysis.m` selects inputs, evaluates cases, records responses, integrates the numerical reference, measures responses, writes CSVs, and exports figures. It sets axes, label, legend, and background colors explicitly. A white figure background alone does not force dark text.

PNG is a raster image. A `.fig` file retains MATLAB graphics objects for editing. CSV stores data, not plots. Markdown holds prose and equations. The optional Python renderer creates the checked-in visual style from the nominal equations; it is documented as such and is not independent physical validation.

Git records file changes in commits. A branch isolates proposed work; a pull request presents the difference; the automated workflow runs a clean checkout in MATLAB. A successful status applies to the tested commit and environment. Its downloadable artifact contains run evidence. After remote updates, fetch and pull before editing locally to avoid accidentally restoring older files.

**Exercise:** from a clean extracted copy, run the tests and analysis, identify each generated output, and explain which files belong in Git and why generated binary `.fig` files are available as artifacts rather than committed by default.

## 12. Questions to practice aloud

Answer each in two layers: a plain-language explanation, then the supporting equation or code reference.

1. What is one degree of freedom, and why do you have two states?
2. Why is the damper force proportional to velocity? What alternative friction law would change the model?
3. Where did gravity go?
4. Why do you need the sine term when initial velocity is zero?
5. What changes mathematically at critical damping?
6. Why can increasing damping increase settling time?
7. What exactly does $\zeta=0.780$ satisfy, and which assumptions does it require?
8. Why does the settling plot have sharp drops?
9. Why does every case have the same initial acceleration?
10. How can displacement grow during part of a cycle while energy decreases?
11. What does a $10^{-11}$ m analytical/numerical discrepancy demonstrate?
12. How could a coarse time grid fool the settling helper?
13. How do the nonzero-initial-velocity tests strengthen the project?
14. How would you identify $m,k,c$ from a real system? Which measurements would be insufficient alone?
15. What would have to change before calling this an aircraft-wing model?
16. Which parts of the workflow did you implement, derive, adapt, debug, or receive help with?

For question 14, a free-decay frequency and decay rate primarily identify parameter ratios, not all three dimensional parameters uniquely. You need an additional known scale, such as mass or a static stiffness measurement. This follows from the scaling invariance in Section 2.

For question 16, describe your actual contribution and assistance accurately. Being able to explain and maintain the resulting code is the goal; do not imply that you independently authored parts you received assistance with.

## 13. Intensive practice sequence

Use approximately 60–90 minutes per session; progress by demonstrated understanding, not hours spent.

| Session | Work | Exit criterion |
| --- | --- | --- |
| 1 | Free-body diagram, assumptions, units, equilibrium shift, nondimensionalization | Derive the ODE and predict parameter scaling without notes. |
| 2 | Roots, all three solutions, initial conditions, energy balance | Derive the release response and energy law on paper. |
| 3 | Reverse peak, 2% design, settling jumps, uncertainty | Reproduce the design ratio and explain why it needs margin. |
| 4 | MATLAB syntax and reusable response function | Write a small analysis script from a blank file and explain every array shape. |
| 5 | State-space, ode45, expm, stable evaluation | Rebuild the numerical comparison and diagnose an intentional sign error. |
| 6 | Sampling, tests, plots, clean-checkout reproduction | Reproduce outputs and explain limits of every verification claim. |
| 7 | Oral defense plus one independent modification | Give a five-minute walkthrough and demonstrate a prediction you verified. |

A good independent modification is an initial-velocity study using `free_response` and `ode45`, with plots and a short explanation of why the original reverse-peak formula no longer applies. Another is comparing 1%, 2%, and 5% settling bands while keeping the system fixed. Neither requires adding Simulink or claiming a controller.

Keep a short engineering notebook: prediction, changed input, observed result, explanation, and any discrepancy. That record demonstrates understanding more convincingly than more software names on a resume.
