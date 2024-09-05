import numpy as np
import main
import matplotlib.pyplot as plt
import random
import tkinter as tk
from tkinter import Canvas
import pandas as pd
from main import s
import math
from tkinter import *
from PIL import Image, ImageTk, EpsImagePlugin


path =  r"D:\\Recrystallization final year project\\Recrys_FYP_2023-24\\energy_misorientation_IQ.txt"  #D:\Python Codes\Recrystallization FYP\HR.xlsx

df = pd.read_csv(path)

df = df.to_numpy()

for i in range(0,len(df[:,0])):
    if df[i+1,0] - df[i,0] != 0:
        stepsize_x = df[i+1,0] - df[i,0]
        break
    
for i in range(0,len(df[:,1])):
    if df[i+1,1] - df[i,1] != 0:
        stepsize_y = df[i+1,1] - df[i,1]
        break
print(stepsize_x,stepsize_y)

df[:,0 ] = (df[:, 0] / (stepsize_x)).astype(int)
df[:,1] = (df[:, 1] / (stepsize_y)).astype(int)

r,c = int(np.max(df[:,0])),int(np.max(df[:,1]))  ##r = x , c= y

init_energy_array = df[:,5]

sigma_energy = np.std(init_energy_array)
mean_energy = np.mean(init_energy_array)
max_energy = np.max(init_energy_array)

print(sigma_energy,mean_energy,max_energy)

nucleation_site = []

for i in df:
    if i[5] >= (mean_energy + 2*sigma_energy):
        nucleation_site.append([i[0],i[1]])
        

number_of_grains = 100

M_m = 10



EA = np.zeros((r+1,c+1,4))
EA = main.s
s = np.zeros((r+1,c+1,3))

for i in df:
    s[int(i[0])][int(i[1])][0] = i[2]  ## theta (average misorientation)
    s[int(i[0])][int(i[1])][1] = i[3]  ## KAM
    s[int(i[0])][int(i[1])][2] = i[5]  ## SE


lattice_angle = np.zeros((r+1,c+1,1))
for i in df:
    lattice_angle[int(i[0])][int(i[1])][0] = i[2]

lattice_status = np.zeros((r+1,c+1),object) ## Stores the color for display makes matrix with 0 but can now store any object

class grain:
    def __init__(self,name,eulerangles,color = "red"):
        self.eulerangles = eulerangles
        self.GB = []
        self.newgrainspx = []
        self.name = name
        self.color = color    
    def isGB(self,x,y):
        #if fetchEA(x,y) == self.eulerangles:
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    if fetchEA(x+i,y+j) != fetchEA(x,y):
                        return True
                    else: pass
            return False        
        #else: pass                

    def updateGB(self):
        for i in self.GB:
            if grain.isGB(self,i[0],i[1]):
                pass
            else: 
                self.GB.remove(i)
        for i in self.newgrainspx:
            #if grain.isGB(self,i[0],i[1]):
            self.GB.append(i)
            #else: pass  
        #self.newgrainspx = []    should be in the final code   


def fetchEA(x,y):
    a = x%(r+1)
    b = y%(c+1)
    return [EA[a,b,0],EA[a,b,1],EA[a,b,2]]

def mobility(misorientation):
    B = 5
    K = 5
    M_m = 10
    return M_m *(1-(math.exp(-1*B*((misorientation/main.theta_m)**K))))


def del_E(EA_M,EA_1,coords_px):  ## EA_M is the euler angles of the grain, EA_1 the eulerangles of the pixel to be changed
    SE_i = 0
    SE_f = 0
    for x in [-1,0,1]:
        for y in [-1,0,1]:
            if (coords_px[0] + x >= 97) or (coords_px[0] + y >= 97): pass
            else:SE_i = SE_i + main.stored_energy(main.theta(np.matmul(main.g(EA_1[0],EA_1[1],EA_1[2]),np.linalg.inv(main.g(fetchEA(coords_px[0]+x,coords_px[1]+y)[0],fetchEA(coords_px[0]+x,coords_px[1]+y)[1],fetchEA(coords_px[0]+x,coords_px[1]+y)[2]))))) ##Might not work
    for x in [-1,0,1]:
        for y in [-1,0,1]:
            if (coords_px[0] + x >= 97) or (coords_px[0] + y >= 97): pass
            else:SE_f = SE_f + main.stored_energy(main.theta(np.matmul(main.g(EA_M[0],EA_M[1],EA_M[2]),np.linalg.inv(main.g(fetchEA(coords_px[0]+x,coords_px[1]+y)[0],fetchEA(coords_px[0]+x,coords_px[1]+y)[1],fetchEA(coords_px[0]+x,coords_px[1]+y)[2]))))) ##Might not work
    
    return SE_f - SE_i


def probability(del_E, misorientation):
    if del_E <= 0 :
        # print('del_E is negative')
        #print(del_E)
        #print((mobility(misorientation)*main.stored_energy(misorientation)*2)/(M_m*main.sigma_m))
        return (mobility(misorientation)*main.stored_energy(misorientation)*2)/(M_m*main.sigma_m) ##function stored energy halves the value
        #return 0.8
    else:
        #print("denominator: ")
        #print(M_m*main.sigma_m)
        return (mobility(misorientation)*main.stored_energy(misorientation)*2)*(np.exp(-1*del_E))/(M_m*main.sigma_m) ## kT term not added
        #return 0.8

def state_change(grain,coords_px):
    pixel_state_initial = fetchEA(coords_px[0],coords_px[1])
    x = probability(del_E(grain.eulerangles, pixel_state_initial,coords_px),np.degrees(main.theta(np.matmul(main.g(grain.eulerangles[0],grain.eulerangles[1],grain.eulerangles[2]),np.linalg.inv(main.g(pixel_state_initial[0],pixel_state_initial[1],pixel_state_initial[2]))))))
    #print(x)
    if random.uniform(0, 1) <= x:
        EA[coords_px[0]%(r+1),coords_px[1]%(c+1),0] = grain.eulerangles[0] # Mod to wrap around
        EA[coords_px[0]%(r+1),coords_px[1]%(c+1),1] = grain.eulerangles[1]
        EA[coords_px[0]%(r+1),coords_px[1]%(c+1),2] = grain.eulerangles[2]
        grain.newgrainspx.append([coords_px[0]%(r+1),coords_px[1]%(c+1)])

        lattice_status[coords_px[0]%(r+1),coords_px[1]%(c+1)] = grain.color  ## Mod to wrap around 
    else: pass    


def print_euler_angles():
    with open(f"sim_output_n={number_of_grains}.txt", "w") as f:
        f.write("phi1,phi,phi2,X,Y,IQ\n")

    for x in range(0, r+1):
        for y in range(0, c+1):
            with open(f"sim_output_n={number_of_grains}.txt", "a") as f:
                f.write("%s,%s,%s,%s,%s,%s\n"%(EA[x, y, 0], EA[x, y, 1], EA[x, y, 2],x*stepsize_x, y*stepsize_y,60))


def generate_random_color():
    # Generate random RGB values
    red = random.randint(0, 255)
    green = random.randint(0, 255)
    blue = random.randint(0, 255)

    # Convert to hexadecimal and format as a color string
    color_string = "#{:02X}{:02X}{:02X}".format(red, green, blue)

    return color_string

def save_canvas_image():
    # Create a PostScript file from the canvas
    canvas.postscript(file="output_image.eps", colormode='color')

    # Use Pillow (PIL) to convert the PostScript file to an image (e.g., PNG)
    EpsImagePlugin.gs_windows_binary = r'C:\Program Files (x86)\gs\gs10.02.1\bin\gswin32c.exe'  # Set the Ghostscript executable path
    img = Image.open("output_image.eps")
    img.save(f"sim_output2_n={number_of_grains}.png", format="png")
    img.close()


# def save_canvas_image():
#     # Get the coordinates and dimensions of the canvas
#     x = root.winfo_rootx() + canvas.winfo_x()
#     y = root.winfo_rooty() + canvas.winfo_y()
#     width = canvas.winfo_width()
#     height = canvas.winfo_height()

#     # Capture the content of the canvas as an image
#     image = ImageGrab.grab(bbox=(x, y, x + width, y + height))

#     # Save the image to a file (e.g., PNG)
#     image.save(f"sim_output_n={number_of_grains}.png", format="png")
#################################################################################################################################################################
# def update_display():
#     canvas.delete("all")
#     min_value = lattice_angle.min()
#     max_value = lattice_angle.max()
    
#     for i in range(r + 1):
#         for j in range(c + 1):
#             value = lattice_angle[i, j]
#             # Map the value to a color gradient
#             normalized_value = (value - min_value) / (max_value - min_value)
#             color = gradient_color(normalized_value)
#             canvas.create_rectangle(i * pixel_size, j * pixel_size, (i + 1) * pixel_size, (j + 1) * pixel_size, fill=color)
#     canvas.update()

def update_display():
    canvas.delete("all")
    for i in range(r):
        for j in range(c):
            if lattice_status[i, j] == 0: 
                color = 'white' 
            else:
                color= lattice_status[i,j]
            canvas.create_rectangle(i * pixel_size, j * pixel_size, (i+1) * pixel_size, (j+1) * pixel_size, fill=color)
    canvas.update()

###############################################################################################################################################################


grains = []




def monte_carlo_step(n=30):
    m=0

    while m < n:
        for i in grains:
            print(i.name)
            for j in i.GB:
                print(j)
                for x in [-1,0,1]:
                    #print(x)
                    for y in [-1,0,1]:
                        #print(y)
                        if lattice_status[(j[0]+x)%(r+1),(j[1]+y)%(c+1)] == 0:
                            if i.eulerangles != fetchEA(j[0]+x,j[1]+y):
                                #print(fetchEA(j[0]+x,j[1]+y))
                                state_change(i,[j[0]+x,j[1]+y])
                                #print(fetchEA(j[0]+x,j[1]+y))
                            
                                #print(i.name)
                                #print(i.GB)
                                #print("GB updated")
            i.updateGB()
        m +=1 
    update_display()


for i in range(1, number_of_grains + 1):
    a= np.random.randint(0,len(nucleation_site))
    nuclii_x = int(nucleation_site[a][0])
    nuclii_y = int(nucleation_site[a][1])
    obj_name = f"grain {i}"
    new_object = grain(obj_name,fetchEA(nuclii_x,nuclii_y),color= generate_random_color())
    new_object.GB.append([nuclii_x,nuclii_y])
    grains.append(new_object)
    lattice_status[nuclii_x,nuclii_y] = new_object.color
# print(grains[0].GB)
# print(grains[0].newgrainspx)
# print(grains[0].eulerangles)
 
#print(grains)
#print(grains[0].GB)
# print(lattice_status[grains[0].GB[0][0],grains[0].GB[0][1]])
# print(lattice_status[grains[0].GB[0][0]+1,grains[0].GB[0][1]+1])

  

root = tk.Tk()
root.title("Random nucleation")

pixel_size = 10
canvas = Canvas(root, width=r*pixel_size, height=c*pixel_size)
canvas.pack()

step_button = tk.Button(root, text="Monte Carlo Step", command=monte_carlo_step)
print_button = tk.Button(root, text= "print Euler angles", command = print_euler_angles)
save_button = tk.Button(root, text= "Save image", command = save_canvas_image)
step_button.pack(side=tk.LEFT, padx=5)
print_button.pack(side=tk.LEFT, padx=5)
save_button.pack(side=tk.LEFT, padx=5)

update_display()

root.mainloop()











#monte_carlo_step()

#


#print(grains[0].isGB(grains[0].newgrainspx[0][0],(grains[0].newgrainspx[0][1])))
#print(fetchEA(grains[0].newgrainspx[0][0],(grains[0].newgrainspx[0][1])))
#print(grains[0].eulerangles)
# def gradient_color(value):
#     # Map a value between 0 and 1 to a color gradient (e.g., from blue to red)
#     # You can customize this function to define your color gradient
#    gray_value = int(255 * (1 - value))
#    return f'#{gray_value:02X}{gray_value:02X}{gray_value:02X}'

# root = tk.Tk()
# root.title("Grain growth simulation")


# pixel_size = 10
# canvas = Canvas(root, width=(r+1)*pixel_size, height=(c+1)*pixel_size)
# canvas.pack()

# update_display()

# root.mainloop()

# if __name__ == "__main__":
#     #get canvas
#     #get random nuclii
#     #grow the nuclii to a certain radii
#     # keep growing the nuclii   