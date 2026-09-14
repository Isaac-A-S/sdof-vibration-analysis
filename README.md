# Transient Response and Damping Design of an SDOF Oscillator


Analytical and numerical study of how viscous damping changes the free response of a linear mass–spring–damper system.


## Overview


This project models the post-disturbance motion of a single-degree-of-freedom (SDOF) oscillator,


$$m\ddot{x}+c\dot{x}+kx=0,$$


and compares underdamped, critically damped, and overdamped behavior. The work combines closed-form solutions, nondimensional analysis, response metrics, and numerical verification with MATLAB.


The model can represent an isolated structural mode in many systems, including vibration isolators, suspensions, machinery, and simplified aerospace structures. The numerical parameters used here are illustrative; this repository does not claim to predict the response of a specific aircraft or passenger comfort.


## Engineering objectives


- Derive the governing equation from a free-body model.


- Express the dynamics using natural frequency and damping ratio.


- Implement exact free-response solutions for all damping regimes.


- Quantify settling time, reverse peak, and energy dissipation.


- Verify the analytical solution against MATLAB's ode45 integrator.


- Show that the preferred damping level depends on the design requirement.


## Model


| Symbol | Quantity | Units |
| --- | --- | --- |
| $m$ | equivalent modal mass | kg |
| $c$ | viscous damping coefficient | N·s/m |
| $k$ | equivalent stiffness | N/m |
| $x(t)$ | displacement from equilibrium | m |
| $x_0$ | initial displacement | m |
| $v_0$ | initial velocity | m/s |


The derived parameters are


$$\omega_n=\sqrt{\frac{k}{m}},\qquad
c_{\mathrm{crit}}=2\sqrt{km},\qquad
\zeta=\frac{c}{c_{\mathrm{crit}}}.$$


The damping ratio $\zeta$ is more useful than $c$ alone because it identifies the response regime independently of the chosen mass and stiffness:


- $0\leq\zeta\<1$: underdamped


- $\zeta=1$: critically damped


- $\zeta>1$: overdamped


For the underdamped case,


$$\omega_d=\omega_n\sqrt{1-\zeta^2}$$


and, for the initial conditions used in this study $(x_0=0.50\ \mathrm{m},\,v_0=0)$,


$$x(t)=x_0e^{-\zeta\omega_nt}
\left[
\cos(\omega_dt)+\frac{\zeta\omega_n}{\omega_d}\sin(\omega_dt)
\right].$$


The critical and overdamped solutions are implemented separately to avoid numerical instability near $\zeta=1$. See [the complete derivation](docs/derivation.md).


## Illustrative parameter set


| Parameter | Value |
| --- | ---: |
| Equivalent mass, $m$ | 1000 kg |
| Equivalent stiffness, $k$ | 40,000 N/m |
| Initial displacement, $x_0$ | 0.50 m |
| Initial velocity, $v_0$ | 0 m/s |
| Natural frequency, $\omega_n$ | 6.3246 rad/s |
| Natural frequency, $f_n$ | 1.0066 Hz |
| Critical damping, $c_{\mathrm{crit}}$ | 12,649.1 N·s/m |


Because these values are not identified from test data, all dimensional results should be interpreted as a controlled case study. The normalized trends in $x/x_0$ versus $\omega_n t$ are the general result.


## Case comparison


The 2% settling time is defined as the earliest time after which


$$|x(t)|\leq 0.02|x_0|$$


remains true.


| $c$ (N·s/m) | $\zeta$ | Regime | $\omega_d$ (rad/s) | 2% settling time (s) | First reverse peak (% of $x_0$) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 2,000 | 0.158 | underdamped | 6.245 | 3.659 | 60.468 |
| 6,000 | 0.474 | underdamped | 5.568 | 1.305 | 18.401 |
| 9,862.6 | 0.780 | underdamped — 2% design | 3.960 | 0.570 | 2.000 |
| 12,000 | 0.949 | underdamped | 2.000 | 0.830 | 0.008 |
| 12,649.1 | 1.000 | critical | — | 0.923 | none |
| 16,000 | 1.265 | overdamped | — | 1.350 | none |


The complete values are stored in [data/case_summary.csv](data/case_summary.csv).


## Key findings


1. **Damping dissipates mechanical energy.** With


$$E=\frac12m\dot{x}^2+\frac12kx^2,$$


the equation of motion gives


$$\frac{dE}{dt}=-c\dot{x}^2\leq0.$$


This provides both a physical interpretation and a validation check for the simulation.


2. **More damping is not always faster.** Increasing $\zeta$ suppresses oscillation, but an overdamped response contains a slow exponential mode. In the selected cases, $c=16{,}000$ N·s/m settles more slowly than the near-critical $c=12{,}000$ N·s/m case.


3. **The design target determines the preferred damping ratio.**


- If sign reversal is prohibited, $\zeta=1$ is the fastest monotonic return to equilibrium.


- If a reverse peak up to 2% is acceptable, the boundary


$$e^{-\pi\zeta/\sqrt{1-\zeta^2}}=0.02$$


gives $\zeta=0.7797$. For the illustrative system, this corresponds to $c=9862.6$ N·s/m and approximately 0.570 s to enter and remain inside the ±2% band.


4. **The initial acceleration exposes the limits of the dimensional example.** Because $v_0=0$,


$$\ddot{x}(0)=-\frac{k}{m}x_0=-20\ \mathrm{m/s^2}.$$


This large value reinforces that the parameters are illustrative rather than calibrated aircraft data.


## Repository structure


~~~text
sdof-vibration-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   └── case_summary.csv
├── docs/
│   └── derivation.md
├── matlab/
│   ├── run_analysis.m
│   ├── free_response.m
│   └── compute_settling_time.m
├── tests/
│   └── run_tests.m
└── figures/README.md
~~~


## Running the analysis


1. Clone or download the repository.


2. Open MATLAB in the repository root.


3. Run:


~~~matlab
run("matlab/run_analysis.m")
~~~


The script calculates all damping cases, compares each analytical response with ode45, writes the case table, and exports publication-ready figures to the figures directory. No add-on toolbox is required.


Run the validation checks with:


~~~matlab
run("tests/run_tests.m")
~~~


## Verification strategy


The repository checks more than whether a curve looks reasonable:


- the analytical solution satisfies the specified initial displacement and velocity;


- analytical displacement agrees with an independent ode45 integration;


- total mechanical energy is non-increasing for $c\geq0$;


- underdamped, critical, and overdamped branches remain continuous near $\zeta=1$;


- every reported settling time uses the same explicit ±2% definition.


## Engineering scope and limitations


This is a linear, lumped-parameter, free-response model. It assumes one translational degree of freedom, constant $m$, $c$, and $k$, small displacement, and viscous damping. A real flexible aircraft includes multiple bending and torsional modes, rigid-body coupling, aerodynamic forces, and time-dependent gust excitation. NASA flexible-aircraft studies use multi-degree-of-freedom or distributed structural models and explicit gust inputs rather than a single initial displacement.


Accordingly, the project demonstrates vibration-analysis methods; it does not establish an aircraft design recommendation or a passenger-comfort prediction.


## Next extensions


- Add a forcing input $F(t)$ for discrete gusts or continuous turbulence.


- Model base excitation rather than initial-displacement-only response.


- Convert the system to state-space form for controls work.


- Identify $m$, $c$, and $k$ from measured free-decay data.


- Extend the model to multiple coupled structural modes.


- Compare passive damping with active feedback control.


- Evaluate acceleration or RMS acceleration only after defining a physically valid output location and comfort criterion.


## References


- E. Kreyszig, *Advanced Engineering Mathematics*, Wiley.


- R. C. Hibbeler, *Engineering Mechanics: Dynamics*, Pearson.


- N. T. Nguyen and I. Tuzcu, [Flight Dynamics of Flexible Aircraft with Aeroelastic and Inertial Force Interactions](https://ntrs.nasa.gov/citations/20100023415), NASA Technical Reports Server.


- C. J. Funk, B. Perry III, and W. A. Silva, [A Summary of Revisions Applied to a Turbulence Response Analysis Method for Flexible Aircraft Configurations](https://ntrs.nasa.gov/citations/20140011901), NASA Technical Reports Server.


- MathWorks, [ode45 documentation](https://www.mathworks.com/help/matlab/ref/ode45.html).


## Author


Isaac — mechanical engineering student focused on dynamics, controls, and aerospace systems.