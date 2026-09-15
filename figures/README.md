# Figure provenance

The three committed PNGs are high-contrast presentation renders from
`tools/render_clean_figures.py` (NumPy/Matplotlib), using the nominal model,
case parameters, and sampled settling definition. No measured data are used.

`matlab/run_analysis.m` independently exports the same three plot types with
explicit dark text, white backgrounds, and corresponding editable `.fig`
files. Its numerical outputs and rendering are exercised by the MATLAB
workflow. Presentation typography may differ between Python and MATLAB.

- `free_response_comparison.png`: six normalized free responses and ±2% band.
- `energy_decay.png`: normalized mechanical energy on a logarithmic scale.
- `settling_time_sweep.png`: finite-grid settling estimate across damping ratios.

The Python renderer is a presentation utility, not a separate validation of
MATLAB. The MATLAB tests use `ode45`, `expm`, and physical/mathematical checks.
