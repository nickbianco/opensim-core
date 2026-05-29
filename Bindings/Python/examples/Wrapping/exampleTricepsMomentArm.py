# -------------------------------------------------------------------------- #
#                  OpenSim:  exampleTricepsMomentArm.py                      #
# -------------------------------------------------------------------------- #
# The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  #
# See http://opensim.stanford.edu and the NOTICE file for more information.  #
# OpenSim is developed at Stanford University and supported by the US        #
# National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    #
# through the Warrior Web program.                                           #
#                                                                            #
# Copyright (c) 2005-2025 Stanford University and the Authors                #
# Author(s): Nicholas Bianco                                                 #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the 'License') you may     #
# not use this file except in compliance with the License. You may obtain a  #
# copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an 'AS IS' BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
# -------------------------------------------------------------------------- #

# This script builds a minimal right-arm model (scapula_r, humerus_r, ulna_r)
# extracted from athlete.osim and reproduces the moment-arm sign flip observed
# in the triceps long head at the elbow joint. Body properties, joint frames,
# wrapping obstacle geometry, and path points all match the full model.

import os
import numpy as np
import matplotlib.pyplot as plt
import opensim as osim

# =========================================================================== #
# Reference: athlete.osim (all joints held at default except elbow_flex_r)
# =========================================================================== #

model = osim.Model(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'athlete.osim'))
state = model.initSystem()

elbow_flex = model.updCoordinateSet().get('elbow_flex_r')
triceps = osim.Millard2012EquilibriumMuscle.safeDownCast(
    model.getMuscles().get('triceps_long_r'))
path = osim.Scholz2015GeometryPath.safeDownCast(triceps.updPath())

elbow_flex.setValue(state, np.radians(30))
model.realizePosition(state)
moment_arm = path.computeMomentArm(state, elbow_flex)
print(f"  elbow_flex_r = {np.degrees(elbow_flex.getValue(state)):.1f} deg: "
      f"computeMomentArm = {moment_arm:+.4f} m")





state = model.initSystem()
elbow_flex.setValue(state, q)
model.realizePosition(state)

# computeCouplingVector
# Zero U, set target coord speed = 1, realize velocity, project constraints.
# For an unconstrained model, coupling = [1] (unit vector for the one DOF).
state.setU(osim.Vector(state.getNU(), 0.0))
elbow_flex.setSpeedValue(state, 1.0)
model.realizeVelocity(state)
# projectU is a no-op here (no kinematic constraints).
coupling = np.array([state.getU().get(i) for i in range(state.getNU())])
coupling /= elbow_flex.getSpeedValue(state)
print(f" coupling = {np.round(coupling, 6)}")

# ── MomentArmSolver line 70: s_ma.updU() = 0 ───────────────────────────
state.setU(osim.Vector(state.getNU(), 0.0))
model.realizeVelocity(state)

# path.addInEquivalentForces(state,1,bF,mF)
# For a PathActuator (optimalForce=1, control=1) in a zero-gravity model,
# getRigidBodyForces() after realizeDynamics() == addInEquivalentForces
# with tension=1 applied to the same state.
model.realizeDynamics(state)
bf = model.getRigidBodyForces(state)

# Mobilized-body-index → body name.
mbod_name = {0: 'Ground'}
for bi in range(model.getBodySet().getSize()):
    body = model.getBodySet().get(bi)
    mbod_name[int(body.getMobilizedBodyIndex())] = body.getName()

# Snapshot body forces before state modification.
# VectorOfSpatialVec: len(bf) gives count; sv.get(0)=torque Vec3,
# sv.get(1)=force Vec3; Vec3.to_numpy() converts to ndarray.
body_F = {}
for idx in range(len(bf)):
    sv = bf.get(idx)
    body_F[idx] = (
        sv.get(0).to_numpy(),   # torque τ_B
        sv.get(1).to_numpy(),   # force  F_B
    )
print(f"  [L78]  body forces from addInEquivalentForces (tension=1):")
for idx, (tau, f) in body_F.items():
    if np.linalg.norm(np.concatenate([tau, f])) < 1e-12:
        continue
    nm = mbod_name.get(idx, f'MBod[{idx}]')
    print(f"           {nm:<14}  τ={np.round(tau,5)}  F={np.round(f,5)}")

# ── MomentArmSolver lines 84-85: matter.multiplyBySystemJacobianTranspose
# genForces[i] = Σ_body ( τ_B·ω_B(i) + F_B·v_B(i) )
# where ω_B(i), v_B(i) are body velocities when U[i]=1, all others=0.
# We compute this for i = elbow_flex_r via virtual work.
s.setU(osim.Vector(s.getNU(), 0.0))
coord.setSpeedValue(s, 1.0)   # U = e_i for elbow_flex_r
m.realizeVelocity(s)

print(f"  [L84]  J^T * bodyForces  (per-body virtual-work contribution):")
gen_force = 0.0
for idx, (tau, f) in body_F.items():
    if np.linalg.norm(np.concatenate([tau, f])) < 1e-12:
        continue
    nm = mbod_name.get(idx, f'MBod[{idx}]')
    if nm == 'Ground':
        omega = v = np.zeros(3)
    else:
        bv    = m.getBodySet().get(nm).getVelocityInGround(s)
        omega = bv.get(0).to_numpy()
        v     = bv.get(1).to_numpy()
    tau_w = float(np.dot(tau, omega))
    f_w   = float(np.dot(f, v))
    gen_force += tau_w + f_w
    print(f"           {nm:<14}  "
            f"v={np.round(v,5)}  ω={np.round(omega,5)}  "
            f"τ·ω={tau_w:+.5f}  F·v={f_w:+.5f}  "
            f"subtotal={tau_w+f_w:+.5f}")

# ── MomentArmSolver line 91: return ~coupling * generalizedForces ────────
# (single-DOF model: gen_forces = [gen_force], coupling = [1])
ma_manual = float(np.dot(coupling, np.full(len(coupling), gen_force)))
ok = abs(ma_manual - ma_api) < 1e-3
print(f"  [L91]  coupling · genForces = {ma_manual:+.6f} m"
        f"  ({'== computeMomentArm ✓' if ok else f'MISMATCH vs computeMomentArm: Δ={ma_manual-ma_api:+.5f}'})")








