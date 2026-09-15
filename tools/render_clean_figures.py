"""Render high-contrast figures for the verified SDOF vibration project.

This is a optional presentation renderer. It evaluates the same closed-form
responses, parameters, time grids, and settling-time definition used by the
MATLAB implementation; it does not alter the engineering results.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DIR = Path(__file__).resolve().parents[1] / "figures"

MASS = 1000.0
STIFFNESS = 40000.0
X0 = 0.50
V0 = 0.0
SETTLING_BAND = 0.02

WN = np.sqrt(STIFFNESS / MASS)
C_CRITICAL = 2.0 * np.sqrt(STIFFNESS * MASS)

ALLOWED_REVERSE_PEAK = 0.02
LOG_PEAK = -np.log(ALLOWED_REVERSE_PEAK)
ZETA_DESIGN = LOG_PEAK / np.sqrt(np.pi**2 + LOG_PEAK**2)
C_DESIGN = ZETA_DESIGN * C_CRITICAL

CASE_NAMES = [
    "Low damping",
    "Moderate damping",
    "2% design",
    "Near critical",
    "Critical",
    "Overdamped",
]

C_VALUES = np.array([2000.0, 6000.0, C_DESIGN, 12000.0, C_CRITICAL, 16000.0])
ZETAS = C_VALUES / C_CRITICAL

# Colorblind-friendly palette with enough contrast on a white background.
COLORS = ["#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#009E73", "#6B4C9A"]


def free_response(t: np.ndarray, damping: float) -> tuple[np.ndarray, np.ndarray]:
    """Return exact displacement and velocity for the free SDOF response."""
    zeta = damping / C_CRITICAL
    alpha = zeta * WN

    if zeta < 1.0 - 1e-12:
        wd = WN * np.sqrt(1.0 - zeta**2)
        a_const = X0
        b_const = (V0 + alpha * X0) / wd
        oscillation = a_const * np.cos(wd * t) + b_const * np.sin(wd * t)
        x = np.exp(-alpha * t) * oscillation
        v = np.exp(-alpha * t) * (
            -alpha * oscillation
            - a_const * wd * np.sin(wd * t)
            + b_const * wd * np.cos(wd * t)
        )
    elif abs(zeta - 1.0) <= 1e-12:
        a_const = X0
        b_const = V0 + WN * X0
        x = (a_const + b_const * t) * np.exp(-WN * t)
        v = (b_const - WN * (a_const + b_const * t)) * np.exp(-WN * t)
    else:
        root_term = WN * np.sqrt(zeta**2 - 1.0)
        r1 = -alpha + root_term
        r2 = -alpha - root_term
        c1 = (V0 - r2 * X0) / (r1 - r2)
        c2 = (r1 * X0 - V0) / (r1 - r2)
        x = c1 * np.exp(r1 * t) + c2 * np.exp(r2 * t)
        v = c1 * r1 * np.exp(r1 * t) + c2 * r2 * np.exp(r2 * t)

    return x, v


def settling_time(t: np.ndarray, x: np.ndarray) -> float:
    """Finite-record estimate of entry into the 2% band, with linear interpolation."""
    distance = np.abs(x) - SETTLING_BAND * abs(X0)
    outside = np.flatnonzero(distance > 0.0)
    if outside.size == 0:
        return float(t[0])
    last = int(outside[-1])
    if last == t.size - 1:
        return float("nan")
    t1, t2 = t[last], t[last + 1]
    y1, y2 = distance[last], distance[last + 1]
    return float(t1 - y1 * (t2 - t1) / (y2 - y1))


def style_axes(ax: plt.Axes) -> None:
    """Apply a consistent high-contrast engineering-report style."""
    dark = "#111827"
    ax.set_facecolor("white")
    ax.tick_params(axis="both", colors=dark, labelsize=10, width=0.9)
    ax.xaxis.label.set_color(dark)
    ax.yaxis.label.set_color(dark)
    ax.title.set_color(dark)
    for spine in ax.spines.values():
        spine.set_color("#374151")
        spine.set_linewidth(0.9)
    ax.grid(True, which="major", color="#D1D5DB", linewidth=0.7, alpha=0.85)
    ax.grid(True, which="minor", color="#E5E7EB", linewidth=0.45, alpha=0.55)
    ax.set_axisbelow(True)


def style_legend(legend: plt.Legend) -> None:
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("#9CA3AF")
    frame.set_alpha(1.0)
    for text in legend.get_texts():
        text.set_color("#111827")


def save_figure(fig: plt.Figure, filename: str) -> None:
    fig.savefig(
        OUTPUT_DIR / filename,
        dpi=300,
        facecolor="white",
        edgecolor="white",
        bbox_inches="tight",
        pad_inches=0.12,
        metadata={"Creator": "SDOF Vibration Analysis portfolio project"},
    )
    plt.close(fig)


def render_free_response(t: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 6.2), facecolor="white")
    ax.axhspan(-SETTLING_BAND, SETTLING_BAND, color="#E5E7EB", alpha=0.65, zorder=0)

    for name, damping, zeta, color in zip(CASE_NAMES, C_VALUES, ZETAS, COLORS):
        x, _ = free_response(t, damping)
        width = 2.35 if name == "2% design" else 1.75
        ax.plot(t, x / X0, color=color, linewidth=width,
                label=rf"{name}: $\zeta={zeta:.3f}$")

    ax.axhline(SETTLING_BAND, color="#4B5563", linestyle="--", linewidth=1.0)
    ax.axhline(-SETTLING_BAND, color="#4B5563", linestyle="--", linewidth=1.0)
    ax.axhline(0.0, color="#6B7280", linestyle=":", linewidth=0.9)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(-0.82, 1.05)
    ax.set_xlabel("Time, $t$ (s)", fontsize=11.5)
    ax.set_ylabel("Normalized displacement, $x/x_0$", fontsize=11.5)
    ax.set_title("Free response across damping regimes", fontsize=14, weight="semibold", pad=11)
    style_axes(ax)
    legend = ax.legend(loc="center left", bbox_to_anchor=(1.015, 0.5), fontsize=9.4,
                       frameon=True, borderpad=0.8)
    style_legend(legend)
    fig.subplots_adjust(right=0.74)
    save_figure(fig, "free_response_comparison.png")


def render_energy(t: np.ndarray) -> None:
    initial_energy = 0.5 * MASS * V0**2 + 0.5 * STIFFNESS * X0**2
    fig, ax = plt.subplots(figsize=(11.2, 6.2), facecolor="white")

    for name, damping, zeta, color in zip(CASE_NAMES, C_VALUES, ZETAS, COLORS):
        x, v = free_response(t, damping)
        energy = 0.5 * MASS * v**2 + 0.5 * STIFFNESS * x**2
        normalized = np.maximum(energy / initial_energy, np.finfo(float).tiny)
        width = 2.35 if name == "2% design" else 1.75
        ax.semilogy(t, normalized, color=color, linewidth=width,
                    label=rf"{name}: $\zeta={zeta:.3f}$")

    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(1e-12, 1.25)
    ax.set_xlabel("Time, $t$ (s)", fontsize=11.5)
    ax.set_ylabel("Normalized mechanical energy, $E/E_0$", fontsize=11.5)
    ax.set_title("Mechanical-energy dissipation", fontsize=14, weight="semibold", pad=11)
    style_axes(ax)
    legend = ax.legend(loc="center left", bbox_to_anchor=(1.015, 0.5), fontsize=9.4,
                       frameon=True, borderpad=0.8)
    style_legend(legend)
    fig.subplots_adjust(right=0.74)
    save_figure(fig, "energy_decay.png")


def render_settling_sweep() -> None:
    zeta_sweep = np.linspace(0.05, 2.00, 391)
    t_sweep = np.arange(0.0, 15.0 + 0.0005, 0.001)
    times = np.empty_like(zeta_sweep)

    for index, zeta in enumerate(zeta_sweep):
        response, _ = free_response(t_sweep, zeta * C_CRITICAL)
        times[index] = settling_time(t_sweep, response)

    design_response, _ = free_response(t_sweep, C_DESIGN)
    design_time = settling_time(t_sweep, design_response)

    fig, ax = plt.subplots(figsize=(10.2, 6.2), facecolor="white")
    ax.plot(zeta_sweep, times, color="#0072B2", linewidth=2.0, label="Numerical sweep")
    ax.scatter([ZETA_DESIGN], [design_time], s=72, color="#D55E00", edgecolor="#111827",
               linewidth=0.9, zorder=5, label="2% reverse-peak boundary")
    ax.axvline(1.0, color="#374151", linestyle="--", linewidth=1.25)
    ax.text(1.025, 13.45, r"Critical damping, $\zeta=1$", color="#111827",
            fontsize=9.5, va="top")
    ax.annotate(
        rf"$\zeta={ZETA_DESIGN:.3f}$,  $t_s={design_time:.3f}$ s",
        xy=(ZETA_DESIGN, design_time),
        xytext=(0.58, 2.35),
        arrowprops={"arrowstyle": "->", "color": "#4B5563", "lw": 1.0},
        fontsize=9.5,
        color="#111827",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#9CA3AF"},
    )
    ax.set_xlim(zeta_sweep[0], zeta_sweep[-1])
    ax.set_ylim(0.0, 14.0)
    ax.set_xlabel(r"Damping ratio, $\zeta$", fontsize=11.5)
    ax.set_ylabel("Two-percent settling time (s)", fontsize=11.5)
    ax.set_title(r"Settling-time tradeoff for $x(0)=x_0$ and $\dot{x}(0)=0$",
                 fontsize=14, weight="semibold", pad=11)
    style_axes(ax)
    legend = ax.legend(loc="upper right", fontsize=9.5, frameon=True)
    style_legend(legend)
    save_figure(fig, "settling_time_sweep.png")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.unicode_minus": True,
        "savefig.transparent": False,
    })

    t = np.arange(0.0, 4.0 + 0.00005, 0.0001)
    render_free_response(t)
    render_energy(t)
    render_settling_sweep()

    design_x, _ = free_response(t, C_DESIGN)
    design_ts = settling_time(t, design_x)
    assert abs(ZETA_DESIGN - 0.7797032674) < 1e-9
    assert abs(C_DESIGN - 9862.5528957) < 1e-6
    assert abs(design_ts - 0.569603) < 5e-6

    print(f"zeta_design={ZETA_DESIGN:.9f}")
    print(f"c_design={C_DESIGN:.6f} N*s/m")
    print(f"design_settling_time={design_ts:.6f} s")
    for output in sorted(OUTPUT_DIR.glob("*.png")):
        print(output)


if __name__ == "__main__":
    main()
