import sys
sys.path.insert(0, '../build/')

import ultrasonics as us

# creating a ultrasonics simulation
sim = us.Simulation()

# creating a transducer at [0, 0, 0] with a director of [0, 1, 0]
t = us.Transducer([0, 0, 0], [0, 1, 0])

# adding transducer to simulation
sim.add_transducer(t)

# defining particle mass and diameter and position
pm = 7.176591426e-7
pd = 0.0042
ppos = [0.0001, 0.006, 0]

# creating a particle
p = us.Particle(ppos, pm, pd)

# printing the Gorkov potential for that particle
print("Gorkov potential:", sim.Gorkov_potential(p))

# defining the number of points, lengths, and origin points for the
# Field we are about to create
N = [200, 200, 1]
L = [0.1, 0.05, 0.0]
origin = [-0.05, 0.005, 0.0]

# calculating Gorkov potential field
gorkov = sim.Gorkov_potential_field(N, L, origin, pm, pd)

# outputing the Field to a vtr file viewable with Paraview
us.dump(gorkov, "basic_simulation.vtr")
