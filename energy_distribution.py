import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path =  r"D:\\Recrystallization final year project\\Recrys_FYP_2023-24\\energy_misorientation_IQ.txt"

df = pd.read_csv(path)

df = df.to_numpy()

energy_values = []

for i in df:
    if i[5] == 0:
        pass
    else:
        energy_values.append(i[5])
plt.hist(energy_values, bins=30, color='skyblue', edgecolor='black')  # Adjust bins as needed
plt.xlabel('Misorientation energy (a.u)')
plt.ylabel('Frequency')
plt.title('Energy Distribution for deformed sample')
plt.savefig("Energy Distribution for deformed sample", dpi=1200)
plt.grid(False)  # Add grid lines
plt.show()       
