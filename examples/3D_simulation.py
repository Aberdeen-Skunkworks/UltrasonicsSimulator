# including build directory
import sys
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

# focusing the transducers at the middle of the simulation
focus_point = [origin[0] + L[0] / 2, origin[1] + L[1] / 2, origin[2] + L[2] / 2]
print("focus point:", focus_point)
sim.focus(focus_point)
    
# defining particle mass and diameter
pm = 7.176591426e-7
pd = 0.0042

# defining the number of points for the field we are about to create
N = [20, 20, 20]

# calculating Gorkov potential field
print("Calculating Gorkov potential field")
gorkov = sim.Gorkov_potential_field(N, L, origin, pm, pd)

# outputing the Field to a vtr file viewable with Paraview
print("outputing Gorkov field to file")
us.dump(gorkov, "3D_simulation.vtr")
