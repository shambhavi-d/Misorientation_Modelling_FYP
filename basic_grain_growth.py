import numpy as np
import matplotlib.pyplot as plt
import random
import tkinter as tk
from tkinter import Canvas

# Define parameters
L = 51  # Size of the lattice
T = 2.0  # Temperature (you can adjust this)
n = 10
states = [0,1,2]
# Initialize the lattice with random spins (0 or 1)
lattice = np.random.choice([0], size=(L, L))
lattice_copy = np.random.choice([0], size=(L, L))
nuclii_x = np.random.randint(0,L-1,(n))
nuclii_y =np.random.randint(0,L-1,(n))

for i in range(0,n):
    lattice[nuclii_x[i],nuclii_y[i]] = random.choice([1,2])
lattice_copy = lattice.copy()

# Create a function to update the lattice for one Monte Carlo step
def monte_carlo_step():
    global lattice , lattice_copy
    for x in range(0,L-1,1):
            for y in range(0,L-1,1):
                if calculate_energy(lattice,x,y)=='state_two'and lattice[x,y]==0:
                    lattice_copy[x, y] = 2
                if calculate_energy(lattice,x,y)=='state_one'and lattice[x,y]==0:
                    lattice_copy[x,y] = 1    
                if calculate_energy(lattice,x,y)=='None' and lattice[x,y]==0:
                    lattice_copy[x,y] = 0
                if calculate_energy(lattice,x,y)=='Random' :
                    lattice_copy[x,y] = random.choice([1,2])
    lattice = lattice_copy.copy()
    update_display()
        

# Define a function to calculate the energy of a given spin and its neighbors
def calculate_energy(lattice, x, y):
    #energy = 0
    state_one = 0
    state_two = 0
    for dx in [-1,0, 1]:
        for dy in [-1,0, 1]:
            if lattice[x + dx, y + dy] == 1:
                state_one +=1
            if lattice[x + dx, y + dy] == 2:
                state_two+=1    
    return 'state_two' if (state_two > state_one) else 'state_one' if (state_one > state_two) else 'None' if ((state_one == 0) and (state_two == 0)) else 'Random'

# Create a function to update the display
def update_display():
    canvas.delete("all")
    for i in range(L):
        for j in range(L):
            color = 'white' if lattice[i, j] == 0 else 'Red' if lattice[i,j] == 1 else 'Green' 
            canvas.create_rectangle(i * pixel_size, j * pixel_size, (i+1) * pixel_size, (j+1) * pixel_size, fill=color)
    canvas.update()

# Create the main tkinter window
root = tk.Tk()
root.title("2D Ising Model")

# Create a canvas for displaying the lattice
pixel_size = 10
canvas = Canvas(root, width=L*pixel_size, height=L*pixel_size)
canvas.pack()

# Create a button to perform a Monte Carlo step
step_button = tk.Button(root, text="Monte Carlo Step", command=monte_carlo_step)
step_button.pack()

# Initialize the display
update_display()

# Start the tkinter main loop
root.mainloop()

# for x in range(0,L,1):
#     for y in range(0,L,1):
#         for dx in [-1, 1]:
#             for dy in [-1, 1]:
#                 print([x+dx,y+dy])