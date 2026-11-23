# rigidFlexibleManipulator

The rigid–flexible manipulator of Ibrahimbegović & Mamouri (2000), *On rigid
components and joint constraints in nonlinear dynamics of flexible multibody
systems employing 3D geometrically exact beam model*, CMAME 188, 805–831,
section 6.2 (Figs. 9–11). It is run with beamFoam alone (`rigidBodyEnd`) and
with moorFV, each with partitioned and monolithic coupling. Results are
compared with each other and with a rigid-link reference.

```sh
./Allrun [endTime]        # beamFoam and moorFV runs side by side, then plots
beamFoam/Allrun [endTime] # beamFoam only
moorFV/Allrun [endTime]   # moorFV only
./plotComparison.py       # plots and table only (plots/); needs numpy and
                          # matplotlib, e.g. ~/venvs/opensees-env/bin/python
./Allclean
```

## The problem (paper, Fig. 9)

- A rigid crank, 0.8 m long, turns about a fixed hinge at 10 rad/s.
- A flexible link, 0.8 m long, is joined to the crank tip by a revolute joint.
  Otherwise the link is free.
- Link properties: EA = GA = 1e6 N, EI = GJ = 1e5 N m², A_ρ = 1 kg/m,
  J_ρ = diag(20, 10, 10) kg m.
- "The length and inertia of each component are the same."
- No gravity, Δt = 0.01 s. The paper plots the free-end displacement over 2 s
  (Fig. 11).

The joint passes no moment, so the link stays straight. Its axial strain is
about 6e-5. The joint follows a circle of radius 0.8 m about the hinge, and
the free end follows it at a nearly constant offset.

## How it is set up here

| Paper | Here |
| --- | --- |
| Hinge, crank tip | Hinge at (−0.8, 0, 0); crank tip (the joint) at the origin; link along +x to (0.8, 0, 0) |
| Crank driven at a constant 10 rad/s | The crank is a free rigid body started at 10 rad/s. Its centre of mass is put at the hinge and held there, and it turns only about z. With Jz = 8.17 kg m² about the hinge (the link's inertia about its end, as "the same inertia" in the paper), the link slows it by under 0.06 % |
| Revolute joint | `joint pin` (force only). The motion is planar, so this is the same as a revolute joint about z |
| Link properties | Set exactly by `0/CQ` (EA, GA), `0/CM` (GJ, EI), `0/ARho` (A_ρ) and `0/CIRho` (J_ρ), which override the circle section in `beamProperties` |
| 4 beam elements | 8 finite-volume cells |
| Start-up not described | The link starts already on the periodic motion; see below |

### Holding the crank on its hinge

- **beamFoam:** `rigidBodyEnd` takes two new optional keys,
  `fixedTranslation yes;` and `rotationAxis (0 0 1);`. Both coupling paths
  already applied these constraints when moorFV supplied them; the keys let a
  beamFoam-only case set them.
- **moorFV:** the six-DoF body uses `point` (centreOfRotation at the hinge) and
  `axis` (0 0 1) constraints.

### The start: why the link is given a small initial spin

The joint force acts at the link's end, not at its centre of mass, so a free
link driven this way cannot keep its direction exactly. With φ the link
direction, a rigid link obeys

```text
φ'' = k sin(ωt − φ),   k = (L/2) m ω² R / (I_c + m L²/4) = 3.133 s⁻²
```

where I_c = A_ρ L³/12 + J_ρ,z L = 8.04 kg m² is the link's inertia about its
centre.

- **Started with no spin, φ'(0) = 0:** the link turns steadily at k/ω =
  0.31 rad/s. By t = 2 s it has turned about 36°, and the free end no longer
  follows the joint. This was checked in a run.
- **Started with φ'(0) = −k/ω = −0.313 rad/s:** the link only wobbles,
  φ ≈ −(k/ω²) sin ωt, about ±0.03 rad. This is the motion in the paper's
  Fig. 11, where the free end follows the joint.

So `0/U`, `0/Accl` and `0/Omega` start the link with the joint velocity
(0, 8, 0) m/s and the spin −0.313 rad/s. The paper's statement that the link
"undergoes a pure translational motion" holds only to within that ±0.03 rad,
and only because J_ρ is large. With the link's own mass inertia alone
(I_c = 0.043 kg m²), the link would visibly swing.

### moorFV

- The six-DoF body runs in `test` mode, so fluid forces are ignored.
- The fluid region is a placeholder: a coarse block round the hinge. Its bottom
  face is the body patch the solver needs.
- `innerDistance` is large, so the whole block turns rigidly with the crank and
  its cells never distort, however far the crank turns.
- Three runs:
  - `partitioned`: `FvBeamNewmark` with 3 PIMPLE outer correctors and
    acceleration relaxation 0.4, as in the other moorFV cases
  - `partitioned_nOuter20`: the same with 20 outer correctors
  - `monolithic`: `beamFoamCoupled`

## Results

Runs of 12 Nov 2026, 2 s at Δt = 0.01 s (`plots/summary.txt`):

| Run | Status | max \|free end − reference\| | max \|free end − beamFoam monolithic\| | Crank speed [rad/s] | Newton iterations per step |
| --- | --- | --- | --- | --- | --- |
| beamFoam, partitioned | finished | 4.4e-3 m | 1.6e-8 m | 9.9947–9.9999 | 16.2 (about 7 coupling iterations) |
| beamFoam, monolithic | finished | 4.4e-3 m | 0 | 9.9947–9.9999 | 4.0 |
| moorFV, partitioned, 3 outer correctors | crashed at 0.67 s | 1.6 m | 1.6 m | 10.4–215 | 11.4 |
| moorFV, partitioned, 20 outer correctors | finished | 1.2e-3 m | 4.5e-3 m | 9.9953–10.0053 | 57.5 |
| moorFV, monolithic | finished | 4.4e-3 m | 1e-14 m | 9.9947–9.9999 | 16.0 (3 outer correctors) |

**Monolithic.**
- beamFoam and moorFV give the same answer to round-off (1e-14 m).
- They need 4 Newton iterations per beam solve throughout. moorFV solves the
  beam once per PIMPLE outer corrector, so it reports 3 × ~5 per step.
- The free end reproduces the paper's Fig. 11: horizontal displacement
  between 0 and −1.6 m, vertical ±0.78 m, period 0.63 s.

**Distance from the reference.**
- The 4.4e-3 m by t = 2 s comes almost entirely from the crank. In these runs
  the crank is a free body, and the link slows it to 9.995 rad/s. The
  reference holds the crank at exactly 10 rad/s, as the paper prescribes.
- The link's direction matches the rigid-link reference throughout
  (`crankAndLink.png`).

**beamFoam partitioned.**
- The coupling iterations converge (about 7 per step, each a beam Newton
  solve).
- The result matches monolithic to 1.6e-8 m, at four times the Newton
  iterations.

**moorFV partitioned, 3 outer correctors** (the usual moorFV setting).
- The body–beam coupling is not converged within each step, and the lag
  speeds the crank up from the first step: 10.40 rad/s after 0.01 s.
- The crank keeps spinning up, to over 200 rad/s, until the beam's Newton
  solve diverges at t = 0.67 s.

**moorFV partitioned, 20 outer correctors.**
- Stable, but still not fully converged: the crank gains about 0.05 % and
  drifts 4.5e-3 m from monolithic by t = 2 s.
- It needs 14 times the Newton iterations of beamFoam monolithic.
- Its smaller distance from the reference is a coincidence: its crank happens
  to stay closer to 10 rad/s.

Run times are not comparable: all five runs ran at once on one machine,
writing every 0.02 s to the same disk.

Figures in `plots/`:
- `tipDisplacement.png`: free end, cf. Fig. 11
- `tipError.png`: free end against the reference and against beamFoam
  monolithic
- `crankAndLink.png`: crank speed and link direction
- `jointForce.png`: joint force, about 64 N
- `iterations.png`: Newton iterations

## Viewing in ParaView

**beamFoam runs** (`beamFoam/partitioned`, `beamFoam/monolithic`):

1. Open `case.foam`, apply **Warp By Vector** with `pointW`, and colour by
   `pointW` or `U` if you like.
2. Open `VTK/rigidBody.vtk.series` in the same session to show the crank: a
   small box at the hinge and an arm to the joint. `beamFoam/Allrun` writes it
   with `tutorials/rigidBodyEndWhirl/writeRigidBodyVTK`. After a partial run,
   rerun it by hand from `beamFoam/`:

   ```sh
   python3 ../../rigidBodyEndWhirl/writeRigidBodyVTK partitioned monolithic --size 0.1
   ```

3. Press Play. Both files are written every 0.02 s, so they animate together.

**moorFV runs** (`moorFV/partitioned`, `partitioned_nOuter20`, `monolithic`):

1. Open `case.foam`. Under **Mesh Regions**, tick `internalMesh` (the fluid
   block, which turns with the crank) and the `beam` region.
2. Apply **Warp By Vector** with `pointW` to the beam block only.

The fluid block's bottom face is the body patch, so the block's rotation shows
the crank angle.

## Files

| Path | Contents |
| --- | --- |
| `beamFoam/base` | beamFoam-only case (`rigidBodyEnd` crank). `beamFoam/Allrun` copies it to `partitioned` and `monolithic` |
| `moorFV/base` | moorFV case (six-DoF crank, beam region `beam`). `moorFV/Allrun` meshes once, then runs `partitioned`, `partitioned_nOuter20` and `monolithic` |
| `plotComparison.py` | Reads all runs and writes `plots/*.png` and `plots/summary.txt` |

Outputs read by the plotting script:

- `postProcessing/0/beamDisplacements_right.dat`: free end
- `beamForcesMoments_left.dat`: joint force
- `beamSolverIterations.dat`: Newton iterations per time step
- the crank's state:
  - beamFoam: `postProcessing/rigidBodyEnd/0/rigidBodyEnd.dat`
  - moorFV: `postProcessing/sixDoF_History/0/sixDoFRigidBodyStateFvBeam.dat`
