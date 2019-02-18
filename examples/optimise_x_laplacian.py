import sys
import numpy as np

# including build directory
sys.path.insert(0, '../build/')

# importing module
import ultrasonics as us

# creating a ultrasonics simulation
sim = us.Simulation()

#defining the lengths and origin point of the simulation
L = [0.1, 0.1, 0.1]
origin = [-0.05, 0.005, -0.05]

# defining the number of transducers in the x and z directions
tx, tz = 10, 10

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
# opt_point = [origin[0] + L[0] / 2, 0.12, origin[2] + L[2] / 2]
print("optimisation point:", opt_point)

# optimising
print("Optimising transducer phases")
sim.optimise_Gorkov_laplacian([opt_point], 2e-6, pm, pd, True, False, False, 50, 1e-0)
#sim.optimise_Gorkov_laplacian([opt_point], 2e-8, pm, pd, True, False, False, 200, 10)

# print("phi:", sim.transducer(0).phi)
us.dump(sim, "optimised_phi_x_laplacian.vtu")

# # defining the number of points for the field we are about to create
# # N = [20, 20, 20]
# N = [100, 100, 1]
# L = [0.1, 0.1, 0.0]
# origin = [-0.05, 0.005, 0.0]

# # calculating Gorkov potential field
# print("Calculating Gorkov potential field")
# gorkov = sim.Gorkov_potential_field(N, L, origin, pm, pd)

# # outputing the Field to a vtr file viewable with Paraview
# print("outputing Gorkov field to file")
# us.dump(gorkov, "optimise_x_laplacian.vtr")

opt_phi = []
for t in range(tx*tz):
    opt_phi.append(sim.transducer(t).phi)

sim.focus(opt_point)

for t in range(tx*tz):
    opt_phi[t] -= sim.transducer(t).phi
    if opt_phi[t] < 0:
        opt_phi[t] += 2 * np.pi
print(opt_phi)

normalised_transducers = []
for t, phi in enumerate(opt_phi):
    normalised_transducers.append(us.Transducer(sim.transducer(t).pos, sim.transducer(t).dir, phi))

us.dump(normalised_transducers, "normalised_transducers.vtu")
