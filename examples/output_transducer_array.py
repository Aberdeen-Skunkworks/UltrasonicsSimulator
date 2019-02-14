# including build directory
import sys
sys.path.insert(0, '../build/')

# importing module
import ultrasonics as us

# creating a ultrasonics simulation
sim = us.Simulation()

#defining the lengths and origin point of the simulation
L = [0.2, 0.1, 0.2]
origin = [-0.1, 0.005, -0.1]

# defining the number of transducers in the x and z directions
tx, tz = 20, 20

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

# focusing the transducers at the middle of the simulation and 12 cm above the transducers
focus_point = [origin[0] + L[0] / 2, 0.12, origin[2] + L[2] / 2]
print("focus point:", focus_point)
sim.focus(focus_point)

us.dump(sim, "output_transducer_array.vtu")
