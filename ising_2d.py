import numpy as np
import matplotlib.pyplot as plt
import random
import tkinter as tk
from tkinter import Canvas

# Define parameters
L = 10  # Size of the lattice
T = 2.0  # Temperature (you can adjust this)

# Initialize the lattice with random spins (0 or 1)
lattice = np.random.choice([0, 1], size=(L, L))

# Create a function to update the lattice for one Monte Carlo step
def monte_carlo_step():
   
        x, y = random.randint(0, L-1), random.randint(0, L-1)
        delta_E = calculate_energy(lattice, x, y) - calculate_energy(lattice, x, y)
        if delta_E <= 0 or random.random() < np.exp(-delta_E / T):
            lattice[x, y] = 1 - lattice[x, y]
        update_display()

# Define a function to calculate the energy of a given spin and its neighbors
def calculate_energy(lattice, x, y):
    energy = 0
    for dx in [-1, 1]:
        for dy in [-1, 1]:
            energy += lattice[(x + dx) % L, (y + dy) % L]
    return -2 * lattice[x, y] * energy

# Create a function to update the display
def update_display():
    canvas.delete("all")
    for i in range(L):
        for j in range(L):
            color = 'black' if lattice[i, j] == 0 else 'white'
            canvas.create_rectangle(i * pixel_size, j * pixel_size, (i+1) * pixel_size, (j+1) * pixel_size, fill=color)
    canvas.update()

# Create the main tkinter window
root = tk.Tk()
root.title("2D Ising Model")

# Create a canvas for displaying the lattice
pixel_size = 50
canvas = Canvas(root, width=L*pixel_size, height=L*pixel_size)
canvas.pack()

# Create a button to perform a Monte Carlo step
step_button = tk.Button(root, text="Monte Carlo Step", command=monte_carlo_step)
step_button.pack()

# Initialize the display
update_display()

# Start the tkinter main loop
root.mainloop()