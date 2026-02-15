import numpy as np
import matplotlib.pyplot as plt

# Load data
cppInputs = np.loadtxt("data/computedInputs.csv", delimiter=",")
cppStates = np.loadtxt("data/states.csv", delimiter=",")
cppOutputs= np.loadtxt("data/outputs.csv", delimiter=",")
cppTrajectory= np.loadtxt("data/trajectory.csv", delimiter=",")

# plot the results
plt.figure(figsize=(6,6))
plt.plot(cppOutputs,linewidth=5, label='Controlled trajectory')
plt.plot(cppTrajectory,'r', linewidth=3, label='Desired trajectory')
plt.xlabel('time steps')
plt.ylabel('Outputs')
plt.legend()
plt.show()


plt.figure(figsize=(6,6))
plt.plot(cppInputs,linewidth=4, label='Computed inputs')
plt.xlabel('time steps')
plt.ylabel('Input')
plt.legend()
plt.show()
