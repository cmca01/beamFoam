#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the rigid-flexible manipulator runs (Ibrahimbegovic & Mamouri 2000,
CMAME 188, section 6.2) with each other and with a rigid-link reference.

Runs (cases that have not been run are skipped):
    beamFoam/partitioned             rigidBodyEnd, partitioned coupling
    beamFoam/monolithic              rigidBodyEnd, monolithic coupling
    moorFV/partitioned               FvBeamNewmark, 3 PIMPLE outer correctors
    moorFV/partitioned_nOuter20      FvBeamNewmark, 20 outer correctors
    moorFV/monolithic                beamFoamCoupled

Reference: the link is nearly rigid (axial strain ~6e-5), so its motion is
that of a rigid link pinned to a joint that moves on a circle of radius R at
the crank speed omega:

    phi'' = k sin(omega t - phi),  k = (L/2) m omega^2 R/(I_c + m L^2/4)

with I_c = A_rho L^3/12 + J_rho,z L about the link centre and
phi'(0) = -k/omega (the cases' drift-free start). Integrated here with RK4 at
constant omega = 10 rad/s, as the paper prescribes. In the runs the crank is a
free body (Jz = 8.17 kg m^2) started at 10 rad/s, so its speed varies slightly.

Figures (written to plots/):
    tipDisplacement.png  free-end horizontal and vertical displacement (the
                         paper's Fig. 11)
    tipError.png         distance of the free end from the reference, and
                         from beamFoam monolithic
    crankAndLink.png     crank angular velocity; link direction
    jointForce.png       magnitude of the joint force
    iterations.png       beam Newton iterations per time step, and totals
                         with run times
and a table (plots/summary.txt, also printed).

Usage: ./plotComparison.py [--show]
Needs numpy and matplotlib (for example ~/venvs/opensees-env/bin/python).
"""

import math
import os
import re
import sys

import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PLOT_DIR = os.path.join(HERE, "plots")

# =============================================================
# PROBLEM DATA (as in the cases: see beamFoam/base/constant/beamProperties)
# =============================================================

OMEGA = 10.0        # crank speed [rad/s]
R = 0.8             # crank length, hinge to joint [m]
L = 0.8             # link length [m]
A_RHO = 1.0         # link mass per length [kg/m]
J_RHO_Z = 10.0      # link rotary inertia per length about z [kg m]
HINGE = np.array([-0.8, 0.0])
JOINT0 = np.array([0.0, 0.0])
TIP0 = np.array([L, 0.0])

M_LINK = A_RHO*L
I_C = A_RHO*L**3/12 + J_RHO_Z*L
K = (L/2)*M_LINK*OMEGA**2*R/(I_C + M_LINK*L**2/4)

# =============================================================
# CASES AND STYLES
# =============================================================

# (directory, label, code, coupling, line style)
CASES = [
    ("beamFoam/partitioned", "beamFoam, partitioned", "beamFoam",
     "partitioned", "-"),
    ("beamFoam/monolithic", "beamFoam, monolithic", "beamFoam",
     "monolithic", "-"),
    ("moorFV/partitioned", "moorFV, partitioned (3 outer correctors)",
     "moorFV", "partitioned", (0, (5, 3))),
    ("moorFV/partitioned_nOuter20", "moorFV, partitioned (20 outer correctors)",
     "moorFV", "partitioned", (0, (1, 1.5))),
    ("moorFV/monolithic", "moorFV, monolithic", "moorFV", "monolithic",
     (0, (5, 3))),
]

COLOURS = {"monolithic": "#2a78d6", "partitioned": "#eb6834"}
# Same coupling, second shade, so the two moorFV partitioned runs differ
SHADES = {"moorFV/partitioned_nOuter20": "#a8431a"}
REFERENCE_COLOUR = "0.15"
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|nan|inf"


# =============================================================
# READERS
# =============================================================

def read_table(path):
    """Numeric rows of a whitespace table; header and comment lines skipped."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts or line.lstrip().startswith("#"):
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    return np.array(rows) if rows else np.empty((0, 0))


def read_sixdof(path):
    """Time, rotation z and omega z from sixDoFRigidBodyStateFvBeam.dat."""
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            vectors = re.findall(r"\(([^)]*)\)", line)
            if len(vectors) < 5:
                continue
            t = float(line.split()[0])
            rot = [float(v) for v in vectors[2].split()]
            omega = [float(v) for v in vectors[4].split()]
            out.append([t, rot[2], omega[2]])
    return np.array(out)


def read_log(path):
    info = {"status": "not run", "runtime": math.nan, "time": 0.0}
    if not os.path.isfile(path):
        return info
    info["status"] = "running"
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("Time = "):
                info["time"] = float(line.split()[2])
            elif line.startswith("ExecutionTime"):
                info["runtime"] = float(line.split()[2])
            elif line.startswith("End"):
                info["status"] = "finished"
            elif "FOAM FATAL" in line or "sigFpe" in line:
                info["status"] = "crashed"
    return info


def read_case(directory, code):
    root = os.path.join(HERE, directory)
    pp = os.path.join(root, "postProcessing")
    tip_path = os.path.join(pp, "0", "beamDisplacements_right.dat")
    if not os.path.isfile(tip_path):
        return None

    case = {}
    tip = read_table(tip_path)
    case["t"] = tip[:, 0]
    case["tip"] = tip[:, 1:3]

    force = read_table(os.path.join(pp, "0", "beamForcesMoments_left.dat"))
    case["force_t"] = force[:, 0]
    case["force"] = np.hypot(force[:, 1], force[:, 2])

    if code == "beamFoam":
        body = read_table(os.path.join(pp, "rigidBodyEnd", "0", "rigidBodyEnd.dat"))
        case["body_t"], case["theta"], case["omega"] = body[:, 0], body[:, 19], body[:, 22]
        case["log"] = read_log(os.path.join(root, "log.beamFoam"))
    else:
        body = read_sixdof(os.path.join(
            pp, "sixDoF_History", "0", "sixDoFRigidBodyStateFvBeam.dat"))
        case["body_t"], case["theta"], case["omega"] = body[:, 0], body[:, 1], body[:, 2]
        case["log"] = read_log(os.path.join(root, "log.interFoam"))

    iters_path = os.path.join(pp, "0", "beamSolverIterations.dat")
    if os.path.isfile(iters_path):
        it = read_table(iters_path)
        case["iter_t"], case["nBeamSolves"], case["nNewton"] = it[:, 0], it[:, 2], it[:, 3]

    # Link direction from the free end and the joint; the joint lies on the
    # crank circle at the crank angle (wrapping of the angle does not matter)
    theta = np.interp(case["t"], case["body_t"], np.unwrap(case["theta"]))
    joint = HINGE + R*np.column_stack([np.cos(theta), np.sin(theta)])
    tip_abs = TIP0 + case["tip"]
    d = tip_abs - joint
    case["phi"] = np.arctan2(d[:, 1], d[:, 0])
    return case


# =============================================================
# REFERENCE: RIGID LINK ON A JOINT MOVING AT CONSTANT OMEGA
# =============================================================

def reference(t_end, dt=1e-4):
    def rhs(t, y):
        return np.array([y[1], K*math.sin(OMEGA*t - y[0])])

    n = int(round(t_end/dt))
    t = np.linspace(0, n*dt, n + 1)
    y = np.zeros((n + 1, 2))
    y[0] = [0.0, -K/OMEGA]
    for i in range(n):
        ti, yi = t[i], y[i]
        k1 = rhs(ti, yi)
        k2 = rhs(ti + dt/2, yi + dt/2*k1)
        k3 = rhs(ti + dt/2, yi + dt/2*k2)
        k4 = rhs(ti + dt, yi + dt*k3)
        y[i + 1] = yi + dt/6*(k1 + 2*k2 + 2*k3 + k4)

    phi, dphi = y[:, 0], y[:, 1]
    ddphi = K*np.sin(OMEGA*t - phi)
    c, s = np.cos(OMEGA*t), np.sin(OMEGA*t)
    joint_disp = np.column_stack([R*(c - 1), R*s])
    tip_disp = joint_disp + L*np.column_stack([np.cos(phi) - 1, np.sin(phi)])

    # Joint force on the link: m times the acceleration of its centre
    e = np.column_stack([np.cos(phi), np.sin(phi)])
    n_ = np.column_stack([-np.sin(phi), np.cos(phi)])
    a_joint = -OMEGA**2*R*np.column_stack([c, s])
    a_centre = a_joint + (L/2)*(ddphi[:, None]*n_ - dphi[:, None]**2*e)
    force = M_LINK*np.hypot(a_centre[:, 0], a_centre[:, 1])
    return {"t": t, "tip": tip_disp, "phi": phi, "force": force}


# =============================================================
# PLOTS
# =============================================================

def style(directory, coupling, ls):
    return dict(color=SHADES.get(directory, COLOURS[coupling]), linestyle=ls,
                linewidth=1.4)


def finish(fig, name):
    fig.tight_layout()
    os.makedirs(PLOT_DIR, exist_ok=True)
    path = os.path.join(PLOT_DIR, name)
    fig.savefig(path, dpi=200)
    print("Saved:", os.path.relpath(path, HERE))


def plot_tip(cases, ref):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    for ax, i, name in zip(axes, [0, 1], ["Horizontal", "Vertical"]):
        ax.plot(ref["t"], ref["tip"][:, i], color=REFERENCE_COLOUR, lw=2.6,
                alpha=0.35, label="Rigid-link reference")
        for (d, label, code, coupling, ls), c in cases:
            ax.plot(c["t"], c["tip"][:, i], label=label, **style(d, coupling, ls))
        ax.set_ylabel(f"{name} free-end displacement [m]")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="lower left")
    axes[1].set_xlabel("Time [s]")
    fig.suptitle("Rigid-flexible manipulator: free-end displacement "
                 "(cf. Ibrahimbegovic & Mamouri 2000, Fig. 11)", fontsize=10)
    finish(fig, "tipDisplacement.png")


def plot_errors(cases, ref, baseline):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    for (d, label, code, coupling, ls), c in cases:
        r = np.column_stack([np.interp(c["t"], ref["t"], ref["tip"][:, i])
                             for i in range(2)])
        axes[0].semilogy(c["t"], np.linalg.norm(c["tip"] - r, axis=1) + 1e-16,
                         label=label, **style(d, coupling, ls))
        if baseline is not None and c is not baseline:
            b = np.column_stack([np.interp(c["t"], baseline["t"], baseline["tip"][:, i])
                                 for i in range(2)])
            axes[1].semilogy(c["t"], np.linalg.norm(c["tip"] - b, axis=1) + 1e-16,
                             label=label, **style(d, coupling, ls))
    axes[0].set_ylabel("|free end - reference| [m]")
    axes[1].set_ylabel("|free end - beamFoam monolithic| [m]")
    axes[1].set_xlabel("Time [s]")
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    finish(fig, "tipError.png")


def plot_crank_link(cases, ref):
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    axes[0].axhline(OMEGA, color=REFERENCE_COLOUR, lw=2.6, alpha=0.35,
                    label=f"Prescribed in the paper ({OMEGA:g} rad/s)")
    axes[1].plot(ref["t"], ref["phi"], color=REFERENCE_COLOUR, lw=2.6, alpha=0.35,
                 label="Rigid-link reference")
    for (d, label, code, coupling, ls), c in cases:
        axes[0].plot(c["body_t"], c["omega"], label=label, **style(d, coupling, ls))
        axes[1].plot(c["t"], c["phi"], label=label, **style(d, coupling, ls))
    # Zoom on the runs that stay near 10 rad/s; name any that leave the band
    off = [(label, c["omega"].max()) for (d, label, *_), c in cases
           if np.max(np.abs(c["omega"] - OMEGA)) > 0.05*OMEGA]
    if off:
        near = [c["omega"] for _, c in cases
                if np.max(np.abs(c["omega"] - OMEGA)) <= 0.05*OMEGA]
        lo = min([OMEGA] + [w.min() for w in near])
        hi = max([OMEGA] + [w.max() for w in near])
        pad = max(0.2*(hi - lo), 1e-3*OMEGA)
        axes[0].set_ylim(lo - pad, hi + pad)
        axes[0].text(
            0.01, 0.04, "Off scale: " + "; ".join(
                f"{label} (up to {w:.0f} rad/s)" for label, w in off),
            transform=axes[0].transAxes, fontsize=8, color="0.3")
    axes[0].set_ylabel("Crank angular velocity [rad/s]")
    axes[1].set_ylabel("Link direction phi [rad]")
    axes[1].set_xlabel("Time [s]")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    finish(fig, "crankAndLink.png")


def plot_force(cases, ref):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(ref["t"], ref["force"], color=REFERENCE_COLOUR, lw=2.6, alpha=0.35,
            label="Rigid-link reference")
    for (d, label, code, coupling, ls), c in cases:
        ax.plot(c["force_t"], c["force"], label=label, **style(d, coupling, ls))
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("|Joint force| [N]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    finish(fig, "jointForce.png")


def plot_iterations(cases):
    have = [(c, case) for c, case in cases if "nNewton" in case]
    if not have:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    for (d, label, code, coupling, ls), c in have:
        axes[0].plot(c["iter_t"], c["nNewton"], label=label, **style(d, coupling, ls))
    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel("Beam Newton iterations per time step")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=7)

    def short(label, c):
        text = (label.replace(", ", "\n", 1).replace("partitioned (", "part. (")
                .replace(" outer correctors)", " outer)").replace(" (", "\n("))
        if c["log"]["status"] == "crashed":
            text += f"\ncrashed at {c['t'][-1]:g} s"
        return text

    labels = [short(lab, c) for (_, lab, *_), c in have]
    totals = [c["nNewton"].sum() for _, c in have]
    colours = [SHADES.get(d, COLOURS[cp]) for (d, _, _, cp, _), _ in have]
    x = np.arange(len(have))
    bars = axes[1].bar(x, totals, color=colours)
    for b, (_, c) in zip(bars, have):
        rt = c["log"]["runtime"]
        axes[1].annotate(f"{rt:.0f} s" if np.isfinite(rt) else "",
                         (b.get_x() + b.get_width()/2, b.get_height()),
                         ha="center", va="bottom", fontsize=8, color="0.3")
    axes[1].set_xticks(x, labels, fontsize=7)
    axes[1].set_ylabel("Total beam Newton iterations (run time above)")
    axes[1].grid(alpha=0.3, axis="y")
    finish(fig, "iterations.png")


# =============================================================
# SUMMARY TABLE
# =============================================================

def summary(cases, ref, baseline):
    head = (f"{'Case':44s} {'status':>9s} {'t end':>6s} {'max|tip-ref|':>13s} "
            f"{'max|tip-mono|':>14s} {'omega range':>17s} {'Newton':>7s} "
            f"{'per step':>9s} {'run time':>9s}")
    lines = [
        "Rigid-flexible manipulator (Ibrahimbegovic & Mamouri 2000, s. 6.2)",
        f"Rigid-link reference: k = {K:.5g} 1/s^2, initial link spin "
        f"{-K/OMEGA:.5g} rad/s",
        "max|tip-mono|: largest distance of the free end from beamFoam "
        "monolithic",
        "", head, "-"*len(head),
    ]
    for (d, label, code, coupling, ls), c in cases:
        r = np.column_stack([np.interp(c["t"], ref["t"], ref["tip"][:, i])
                             for i in range(2)])
        err_ref = np.max(np.linalg.norm(c["tip"] - r, axis=1))
        if baseline is not None:
            b = np.column_stack([np.interp(c["t"], baseline["t"], baseline["tip"][:, i])
                                 for i in range(2)])
            err_mono = np.max(np.linalg.norm(c["tip"] - b, axis=1))
        else:
            err_mono = math.nan
        newton = c["nNewton"].sum() if "nNewton" in c else math.nan
        per_step = c["nNewton"].mean() if "nNewton" in c else math.nan
        lines.append(
            f"{label:44s} {c['log']['status']:>9s} {c['t'][-1]:6.3g} "
            f"{err_ref:13.3e} {err_mono:14.3e} "
            f"{c['omega'].min():8.4f}-{c['omega'].max():<8.4f} "
            f"{newton:7.0f} {per_step:9.2f} {c['log']['runtime']:8.1f}s")
    text = "\n".join(lines) + "\n"
    print(text)
    os.makedirs(PLOT_DIR, exist_ok=True)
    with open(os.path.join(PLOT_DIR, "summary.txt"), "w") as f:
        f.write(text)


def main():
    cases = []
    for spec in CASES:
        case = read_case(spec[0], spec[2])
        if case is None:
            print("Not run, skipped:", spec[0])
            continue
        cases.append((spec, case))
    if not cases:
        sys.exit("No runs found: run ./Allrun first")

    t_end = max(c["t"][-1] for _, c in cases)
    ref = reference(t_end)
    baseline = next((c for s, c in cases if s[0] == "beamFoam/monolithic"), None)

    plot_tip(cases, ref)
    plot_errors(cases, ref, baseline)
    plot_crank_link(cases, ref)
    plot_force(cases, ref)
    plot_iterations(cases)
    summary(cases, ref, baseline)

    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
