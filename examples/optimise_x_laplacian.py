import sys
import numpy as np
import pickle

# including build directory
sys.path.insert(0, '../build/')

# importing module
import ultrasonics as us

#defining the lengths and origin point of the simulation
L = [0.1, 0.1, 0.1]
origin = [-0.05, 0.005, -0.05]

# defining the number of transducers in the x and z directions
tx, tz = 10, 10

# creating a ultrasonics simulation
sim = us.Simulation()

for i in range(tx):
    for j in range(tz):

        # transducer position
        tpos = [origin[0] + L[0] * (i+0.5) / tx,
                0,
                origin[2] + L[2] * (j+0.5) / tz]

        # creating a transducer
        t = us.Transducer(tpos, [0, 1, 0])

        # adding transducer to simulation
        sim.add_transducer(t)

# defining particle mass and diameter
pm = 7.176591426e-7
pd = 0.0042

# optimising the transducer phases for the laplacian in the x
# direction at a point in space
opt_point = [origin[0] + L[0] / 2, origin[1] + L[1] / 2, origin[2] + L[2] / 2]
print("optimisation point:", opt_point)

# optimising
print("Optimising transducer phases")
sim.optimise_Gorkov_laplacian([opt_point], 2e-6, pm, pd, True, False, False, 100, 200, 3e-5, "IPOPT")
# There are a number of possible algorithm labels to use including: IPOPT,
# LFBGS, LN_BOBYQA, and GN_ESCH

us.dump(sim, "optimised_phi_x_laplacian.vtu")

# defining the number of points for the field we are about to create
N = [100, 100, 1]
L = [0.1, 0.1, 0.0]
origin = [-0.05, 0.005, 0.0]

# calculating Gorkov potential field
print("Calculating Gorkov potential field")
gorkov = sim.Gorkov_potential_field(N, L, origin, pm, pd)

# outputing the Field to a vtr file viewable with Paraview
print("outputing Gorkov field to file")
us.dump(gorkov, "optimise_x_laplacian.vtr")

# recreating results of the nature communications paper buy
# normalising the transducer phases by the phases focussed at the
# optimisation point
opt_phi = []
for t in range(tx*tz):
    opt_phi.append(sim.transducer(t).phi)

sim.focus(opt_point)

for t in range(tx*tz):
    opt_phi[t] -= sim.transducer(t).phi
    opt_phi[t] = opt_phi[t] % (2 * np.pi)

normalised_transducers = []
for t, phi in enumerate(opt_phi):
    normalised_transducers.append(us.Transducer(sim.transducer(t).pos, sim.transducer(t).dir, phi))

# outputing the normalised phases to a vtu file
print("outputing normalised phases")
us.dump(normalised_transducers, "normalised_transducers_x.vtu")

