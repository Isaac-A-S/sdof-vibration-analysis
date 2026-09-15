# Mathematical Model and Derivation

## 1. Scope

The system is modeled as a linear single-degree-of-freedom oscillator undergoing free vibration after an initial disturbance. The coordinate $x(t)$ measures displacement from static equilibrium, so constant loads such as weight do not appear in the dynamic equation.

Assumptions:

- one equivalent translational degree of freedom;

- constant equivalent mass $m$, damping coefficient $c$, and stiffness $k$;

- linear elastic restoring force;

- linear viscous damping;

- small motion;

- no external force after $t=0$.

The last assumption is important: this is a post-disturbance free-decay model, not a model of continuous turbulence.

## 2. Equation of motion

The spring and damper oppose displacement and velocity:

$$F_k=-kx,\qquad F_c=-c\dot{x}.$$

Newton's second law gives

$$m\ddot{x}=-c\dot{x}-kx,$$

or

$$\boxed{m\ddot{x}+c\dot{x}+kx=0}.$$

Initial conditions are written generally as

$$x(0)=x_0,\qquad \dot{x}(0)=v_0.$$

## 3. Normalized parameters

Define

$$\omega_n=\sqrt{\frac{k}{m}},\qquad
c_{\mathrm{crit}}=2\sqrt{km},\qquad
\zeta=\frac{c}{c_{\mathrm{crit}}}.$$

Because $c/m=2\zeta\omega_n$, division by $m$ gives

$$\boxed{\ddot{x}+2\zeta\omega_n\dot{x}+\omega_n^2x=0}.$$

Using normalized displacement $y=x/x_0$ and normalized time $\tau=\omega_n t$ gives

$$\frac{d^2y}{d\tau^2}+2\zeta\frac{dy}{d\tau}+y=0.$$

The normalized initial conditions are $y(0)=1$ and $y'(0)=v_0/(\omega_n x_0)$, with the prime taken with respect to $\tau$. For $x_0\ne0$ and fixed normalized initial velocity, the response shape depends only on $\zeta$. The dimensional values of $m$ and $k$ set the time and force scales.

## 4. Characteristic equation

Assume $x=e^{rt}$. Substitution gives

$$mr^2+cr+k=0,$$

or

$$r^2+2\zeta\omega_nr+\omega_n^2=0.$$

Therefore,

$$r_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.$$

The root type produces three response regimes.

## 5. Underdamped response: $0\leq\zeta<1$

Define the damped natural frequency

$$\omega_d=\omega_n\sqrt{1-\zeta^2}.$$

The response is

$$x(t)=e^{-\zeta\omega_nt}
\left[A\cos(\omega_dt)+B\sin(\omega_dt)\right].$$

Applying $x(0)=x_0$ and $\dot{x}(0)=v_0$ gives

$$A=x_0,\qquad
B=\frac{v_0+\zeta\omega_nx_0}{\omega_d}.$$

Hence,

$$\boxed{
x(t)=e^{-\zeta\omega_nt}
\left[
x_0\cos(\omega_dt)+
\frac{v_0+\zeta\omega_nx_0}{\omega_d}\sin(\omega_dt)
\right]
}.$$

For $v_0=0$,

$$x(t)=x_0e^{-\zeta\omega_nt}
\left[
\cos(\omega_dt)+
\frac{\zeta\omega_n}{\omega_d}\sin(\omega_dt)
\right].$$

The exponential factor sets the decay rate. For $0<\zeta<1$, its half-life is

$$t_{1/2}=\frac{\ln 2}{\zeta\omega_n}=\frac{2m\ln 2}{c},$$

but this is the half-life of the exponential envelope, not necessarily the first time the actual displacement reaches $x_0/2$.

## 6. Critically damped response: $\zeta=1$

The characteristic equation has the repeated root $r=-\omega_n$:

$$x(t)=(A+Bt)e^{-\omega_nt}.$$

The initial conditions give

$$A=x_0,\qquad B=v_0+\omega_nx_0.$$

Thus,

$$\boxed{x(t)=
\left[x_0+(v_0+\omega_nx_0)t\right]e^{-\omega_nt}}.$$

When $v_0=0$, this response returns to equilibrium without crossing it. For fixed $m$, $k$, $x_0>0$ and $v_0=0$, it is the fastest monotonic return among the constant nonnegative damping choices. This statement does not apply to arbitrary initial velocity.

## 7. Overdamped response: $\zeta>1$

The roots are real and distinct:

$$r_1=-\zeta\omega_n+\omega_n\sqrt{\zeta^2-1},$$

$$r_2=-\zeta\omega_n-\omega_n\sqrt{\zeta^2-1}.$$

The response is

$$x(t)=C_1e^{r_1t}+C_2e^{r_2t},$$

where

$$C_1=\frac{v_0-r_2x_0}{r_1-r_2},\qquad
C_2=\frac{r_1x_0-v_0}{r_1-r_2}.$$

For the illustrative case $m=1000$ kg, $k=40000$ N/m, $c=16000$ N·s/m, $x_0=0.50$ m, and $v_0=0$,

$$r_1=-3.1010\ \mathrm{s^{-1}},\qquad
r_2=-12.8990\ \mathrm{s^{-1}},$$

and

$$\boxed{x(t)=0.65825e^{-3.1010t}-0.15825e^{-12.8990t}}.$$

The slow $r_1$ mode dominates the long-time response. This is why excessive damping can increase settling time.

## 8. Energy balance

Define the total mechanical energy

$$E(t)=\frac12m\dot{x}^2+\frac12kx^2.$$

Differentiate:

$$\frac{dE}{dt}=m\dot{x}\ddot{x}+kx\dot{x}
=\dot{x}(m\ddot{x}+kx).$$

Using the equation of motion,

$$m\ddot{x}+kx=-c\dot{x},$$

so

$$\boxed{\frac{dE}{dt}=-c\dot{x}^2\leq0}.$$

For positive damping, the damper can only remove mechanical energy. The MATLAB tests use this result as a physical consistency check.

## 9. Short-time series check

The differential equation gives

$$\ddot{x}(0)=-\frac{c}{m}v_0-\frac{k}{m}x_0.$$

For the selected condition $v_0=0$,

$$\ddot{x}(0)=-\omega_n^2x_0.$$

Differentiating the equation once gives

$$x^{(3)}=-\frac{c}{m}\ddot{x}-\frac{k}{m}\dot{x},$$

so

$$x^{(3)}(0)=2\zeta\omega_n^3x_0.$$

The local expansion is therefore

$$\boxed{
x(t)=x_0\left[
1-\frac{\omega_n^2t^2}{2}
+\frac{\zeta\omega_n^3t^3}{3}
+O(t^4)
\right]
}.$$

There is no linear term because $v_0=0$. Damping first appears in the cubic term; it cannot change the instantaneous initial acceleration when the initial velocity is zero.

## 10. Response metrics and damping design

### 10.1 Two-percent settling time

The numerical study defines

$$t_s=\inf\left\{t:\ |x(\tau)|\leq0.02|x_0|
\text{ for every }\tau\geq t\right\}.$$

This explicit definition prevents an early zero crossing from being mistaken for settling.

### 10.2 First reverse peak

For the underdamped response with $v_0=0$, extrema after $t=0$ occur when

$$\sin(\omega_dt)=0.$$

The first reverse peak occurs at

$$t_p=\frac{\pi}{\omega_d}.$$

Its magnitude relative to the initial displacement is

$$\boxed{
M_p=\exp\left(
-\frac{\pi\zeta}{\sqrt{1-\zeta^2}}
\right)
}.$$

If the allowed reverse peak is $\delta$, solve $M_p=\delta$:

$$\boxed{
\zeta_\delta=
\frac{-\ln\delta}
{\sqrt{\pi^2+(\ln\delta)^2}}
}.$$

For $\delta=0.02$,

$$\zeta_{2\%}=0.7797.$$

With $m=1000$ kg and $k=40000$ N/m,

$$c=\zeta c_{\mathrm{crit}}
=0.7797(12649.1)
=9862.6\ \mathrm{N\,s/m}.$$

This design enters the ±2% band in approximately $0.570$ s and reaches a first reverse peak of exactly $-2\%$. By contrast, $\zeta=1$ is preferred when any sign reversal is unacceptable. The result is objective-dependent rather than a universal instruction to maximize damping.

## 11. Illustrative numerical case

For

$$m=1000\ \mathrm{kg},\quad
k=40000\ \mathrm{N/m},\quad
x_0=0.50\ \mathrm{m},\quad
v_0=0,$$

the natural frequency and critical damping are

$$\omega_n=6.3246\ \mathrm{rad/s},$$

$$f_n=\frac{\omega_n}{2\pi}=1.0066\ \mathrm{Hz},$$

$$c_{\mathrm{crit}}=12649.1\ \mathrm{N\,s/m}.$$

The initial acceleration is

$$\ddot{x}(0)=-\omega_n^2x_0=-20\ \mathrm{m/s^2}.$$

This scale is intentionally labeled illustrative. A credible aircraft-specific study would require parameters identified from structural or flight-test data.

## 12. Model limitations

The model excludes:

- distributed mass and stiffness;

- multiple bending and torsional modes;

- rigid-body and aeroelastic coupling;

- aerodynamic damping and stiffness;

- time-varying or stochastic gust forcing;

- nonlinear stiffness, damping, and geometry;

- structural-to-cabin transfer paths;

- a validated human-comfort metric.

The correct interpretation is a study of second-order vibration behavior. A future aircraft application should begin with

$$m\ddot{x}+c\dot{x}+kx=F(t)$$

or a multi-degree-of-freedom/state-space model with physically identified parameters and excitation.

## 13. Numerical evaluation and verification limits

The underdamped implementation evaluates $\sin(\omega_d t)/\omega_d$ as
$t\,\mathrm{sinc}_u(\omega_d t)$, where $\mathrm{sinc}_u(u)=\sin(u)/u$ and
$\mathrm{sinc}_u(0)=1$. The subscript distinguishes it from MATLAB's normalized
`sinc`, which uses $\pi u$. A local series avoids division by a tiny argument.
For overdamping the slow root is rationalized:

$$r_s=-\frac{\omega_n}{\zeta+\sqrt{\zeta^2-1}}.$$

With $g=r_s-r_f>0$ and $B=v_0-r_sx_0$, an equivalent solution is

$$D=e^{r_st}\frac{-\operatorname{expm1}(-gt)}{g},\quad
x=x_0e^{r_st}+BD,\quad v=r_sx+Be^{r_ft}.$$

`expm1(u)` evaluates $e^u-1$ accurately near zero. These expressions reduce
cancellation as the two real roots approach each other. Only an eight-ULP
neighborhood of $\zeta=1$ uses the repeated-root limit; this is a documented
floating-point approximation, not a new physical regime.

The settling helper estimates the last entry from a finite sampled record.
It cannot prove the infinite-time definition by itself. Samples can miss a
narrow excursion. We compare refined grids and longer records for the case
set. For underdamping and $v_0=0$, the decreasing envelope

$$\frac{|x(t)|}{|x_0|}\leq
\frac{e^{-\zeta\omega_n t}}{\sqrt{1-\zeta^2}}$$

bounds the remaining tail. For $\zeta\geq1$, the selected zero-velocity
response is monotone. The helper permits 32 floating-point spacings of the
reference displacement when classifying boundary equality; that allowance
is negligible numerically and is not an engineering safety margin.

At the exact 2% design the first negative peak touches the boundary without
leaving it. A slightly smaller damping ratio makes that peak exceed the band
and delays permanent entry to after the negative peak. Settling time can
therefore jump even though displacement varies smoothly with damping.
A design sitting exactly on this boundary has no allowance for parameter
uncertainty. This repository reports a nominal boundary design, not a robust
optimum for hardware. The formula assumes $v_0=0$ and fixed linear parameters.

Numerical agreement with `ode45` or `expm` verifies implementation of the same
idealized equation. It does not validate the equation against measurements.
The returned acceleration is calculated from the ODE, so testing that same
identity is not independent verification. Central differences of displacement
and velocity provide a separate differential-consistency check.
