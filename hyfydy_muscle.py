import numpy as np
import matplotlib.pyplot as plt

# muscle parameters
penn0 = 0.0
Fmax = 4000.0
lopt = 0.55
lts = 0.25
beta = 0.1

# muscle curve parameters
cT1 = 260.972
cT2 = 7.9706

cL1 = 1.5
cL2 = -2.75
r1 = 0.46899
r2 = 1.80528

cP1 = 1.08027
cP2 = 1.27368

cV1 = 0.227
cV2 = 0.110
Fvmax = 1.6

def tendon_force_length_curve(norm_tendon_length):
    if norm_tendon_length > 1.0:
        tendon_strain = norm_tendon_length - 1.0
        tendon_strain_squared = tendon_strain * tendon_strain
        return cT1 * tendon_strain_squared + cT2 * tendon_strain
    else:
        return 0.0


def active_force_length_curve(norm_fiber_length):
    if (norm_fiber_length > r1) and (norm_fiber_length < r2):
        fiber_strain = norm_fiber_length - 1.0
        fiber_strain_squared = fiber_strain * fiber_strain
        fiber_strain_cubed = fiber_strain_squared * fiber_strain
        return cL1 * fiber_strain_cubed + cL2 * fiber_strain_squared + 1.0
    else:
        return 0.0
    

def passive_force_length_curve(norm_fiber_length):
    if norm_fiber_length > 1.0:
        fiber_strain = norm_fiber_length - 1.0
        fiber_strain_squared = fiber_strain * fiber_strain
        fiber_strain_cubed = fiber_strain_squared * fiber_strain
        return cP1 * fiber_strain_cubed + cP2 * fiber_strain_squared
    else:
        return 0.0


def force_velocity_curve(norm_fiber_velocity):
    if norm_fiber_velocity <= -1:
        return 0.0
    elif norm_fiber_velocity >= 0.0:
        return (Fvmax*norm_fiber_velocity + cV2) / (cV2 + norm_fiber_velocity)
    else:
        return (cV1*(norm_fiber_velocity + 1.0)) / (cV1 - norm_fiber_velocity)
    

# norm_tendon_lengths = np.linspace(0.5, 2.0, 100)
# norm_fiber_lengths = np.linspace(0.5, 2.0, 100)
# norm_fiber_velocities = np.linspace(-1.5, 1.5, 100)

# tendon_forces = [tendon_force_length_curve(l) for l in norm_tendon_lengths]
# active_forces = [active_force_length_curve(l) for l in norm_fiber_lengths]
# passive_forces = [passive_force_length_curve(l) for l in norm_fiber_lengths]
# force_velocities = [force_velocity_curve(v) for v in norm_fiber_velocities]

# fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# axs[0, 0].plot(norm_tendon_lengths, tendon_forces)
# axs[0, 0].set_title('Tendon Force-Length Curve')
# axs[0, 0].set_xlabel('Normalized Tendon Length')
# axs[0, 0].set_ylabel('Tendon Force')

# axs[0, 1].plot(norm_fiber_lengths, active_forces)
# axs[0, 1].set_title('Active Force-Length Curve')
# axs[0, 1].set_xlabel('Normalized Fiber Length')
# axs[0, 1].set_ylabel('Active Force')

# axs[1, 0].plot(norm_fiber_lengths, passive_forces)
# axs[1, 0].set_title('Passive Force-Length Curve')
# axs[1, 0].set_xlabel('Normalized Fiber Length')
# axs[1, 0].set_ylabel('Passive Force')

# axs[1, 1].plot(norm_fiber_velocities, force_velocities)
# axs[1, 1].set_title('Force-Velocity Curve')
# axs[1, 1].set_xlabel('Normalized Fiber Velocity')
# axs[1, 1].set_ylabel('Force Velocity')

# plt.tight_layout()
# plt.show()



def calc_fiber_velocity_solutions(muscle_tendon_length, activation, norm_fiber_length):
    fiber_length = norm_fiber_length * lopt
    tendon_length = muscle_tendon_length - np.sqrt(fiber_length**2 - (lopt * np.sin(penn0))**2)
    norm_tendon_length = tendon_length / lts
    cosPenn = (muscle_tendon_length - tendon_length) / fiber_length

    ft = tendon_force_length_curve(norm_tendon_length)
    fl = active_force_length_curve(norm_fiber_length)
    fp = passive_force_length_curve(norm_fiber_length)

    f_target = (ft / cosPenn) - fp

    # norm_fiber_velocity <= -1
    sol1 = f_target / beta
    if sol1 > -1:
        sol1 = np.nan

    # -1 < norm_fiber_velocity < 0
    a2 = -beta
    a1 = activation*fl*cV1 + f_target + beta*cV1
    a0 = activation*fl*cV1 - f_target*cV1

    discriminant = a1**2 - 4*a2*a0
    sqrt_discriminant = np.sqrt(discriminant)
    sol2 = (-a1 + sqrt_discriminant) / (2 * a2)
    sol3 = (-a1 - sqrt_discriminant) / (2 * a2)
    if sol2 <= -1 or sol2 >= 0:
        sol2 = np.nan
    if sol3 <= -1 or sol3 >= 0:
        sol3 = np.nan

    # norm_fiber_velocity >= 0
    a2 = beta
    a1 = activation*fl*Fvmax + beta*cV2 - f_target
    a0 = activation*fl*cV2 - f_target*cV2
    discriminant = a1**2 - 4*a2*a0
    sqrt_discriminant = np.sqrt(discriminant)
    sol4 = (-a1 + sqrt_discriminant) / (2 * a2)
    sol5 = (-a1 - sqrt_discriminant) / (2 * a2)
    if sol4 < 0:
        sol4 = np.nan
    if sol5 < 0:
        sol5 = np.nan

    return np.array([sol1, sol2, sol3, sol4, sol5])

activation = 0.5
muscle_tendon_length = lopt + lts
norm_fiber_lengths = np.linspace(0.75, 1.25, 100)

solutions = np.zeros((len(norm_fiber_lengths), 5))
for i, norm_fiber_length in enumerate(norm_fiber_lengths):
    solutions[i, :] = calc_fiber_velocity_solutions(muscle_tendon_length, activation, norm_fiber_length)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(norm_fiber_lengths, solutions[:, 0], label='norm_fiber_velocity <= -1')
plt.plot(norm_fiber_lengths, solutions[:, 1], label='-1 < norm_fiber_velocity < 0 (sol2)')
plt.plot(norm_fiber_lengths, solutions[:, 2], label='-1 < norm_fiber_velocity < 0 (sol3)')
plt.plot(norm_fiber_lengths, solutions[:, 3], label='norm_fiber_velocity >= 0 (sol4)')
plt.plot(norm_fiber_lengths, solutions[:, 4], label='norm_fiber_velocity >= 0 (sol5)')
plt.title('Fiber Velocity Solutions vs. Normalized Fiber Length')
plt.xlabel('Normalized Fiber Length')
plt.ylabel('Fiber Velocity Solutions')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
plt.axvline(1.0, color='red', linestyle='--', linewidth=0.5, label='lopt')
plt.legend()
plt.grid()
plt.show()
