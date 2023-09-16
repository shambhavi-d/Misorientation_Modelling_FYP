import numpy as np
import pandas as pd


## The first column of the data file will be read as Phi_1, 2nd- Phi, 3rd phi_2, 4th - x, 5th - y, 6and 7th as CI (confidance index)
path =  r"D:\\Recrystallization final year project\\Recrys_FYP_2023-24\\HR_complete.xlsx"  #D:\Python Codes\Recrystallization FYP\HR.xlsx
#print(path)
df = pd.read_excel(path)
# print(df)
df = df.to_numpy()    ##We use numpy arrays since they are easier to manipulate

stepsize = df[1,3] - df[0,3]

df[:,3] = (df[:, 3] / (stepsize)).astype(int)
df[:,4] = (df[:, 4] / (stepsize)).astype(int)

r,c = int(np.max(df[:,3])),int(np.max(df[:,4])) ##Dimensions of the datafile for furthur use

######################## GLOBAL VARIABLES ###############################

theta_m = 15 # 15 degree is the critical value for misorientation

s = np.zeros((r+1,c+1,4))  ## np.zeros makes a nested list (which can be thought of a a matrix of r+1 x c+1 storing 3 values). 
                           ## We are going to store the angles in this such that if i call S[x][y][0] we'll get phi_1 for x,y co-ordinates
                           ## simlarly S[x][y][1] = phi for x,y and S[x][y][2] = phi_2 for x,y

for i in df:
    # print(i)
    s[int(i[3])][int(i[4])][0] = i[0]
    s[int(i[3])][int(i[4])][1] = i[1]
    s[int(i[3])][int(i[4])][2] = i[2]
    s[int(i[3])][int(i[4])][3] = i[6]

G = np.zeros((r + 1, c + 1, 2, 3, 3)) ## Like above in a r+1 x c+1 matrix now we are storing g and g^-1. We do this after caculating g and g^-1

####################### FUNCTIONS #######################################

def g(phi1,phi,phi2):                                          ## Gets the g value for a given phi_1,phi,phi_2
    g_one = np.array([ [np.cos(phi1), np.sin(phi1), 0] ,[-np.sin(phi1), np.cos(phi1), 0] ,[0,0,1]])
    g_two = np.array([ [np.cos(phi2), np.sin(phi2), 0] ,[-np.sin(phi2), np.cos(phi2), 0] ,[0,0,1]])
    g = np.array([ [1,0,0], [0, np.cos(phi), np.sin(phi)] ,[0, -np.sin(phi), np.cos(phi)]])

    b = np.matmul(g,g_one)
    c = np.matmul(g_two,b)

    return c

def theta(del_g): ## Takes del_g (a 3x3 numpy array) and gives back the minimum theta 
    symmetry_matrix = np.array([
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
    [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
    [[0, 1, 0], [0, 0, -1], [-1, 0, 0]],
    [[0, 0, -1], [1, 0, 0], [0, -1, 0]],
    [[0, 0, -1], [-1, 0, 0], [0, 1, 0]],
    [[0, 0, 1], [-1, 0, 0], [0, -1, 0]],
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]]])
    
    min_value = np.inf

    for i in symmetry_matrix:
        n_val = (np.arccos((np.trace(np.matmul( i, del_g) ) - 1) / 2))
        #print(n_val)
        if  n_val < min_value:
            min_value = n_val
     
    return min_value


###########Setting up G ##########################

for x in range(0,r+1):
    for y in range(0,c+1):
        G[x,y]= [g(s[x,y,0] , s[x,y,1] , s[x,y,0]), np.linalg.inv(g(s[x,y,0] , s[x,y,1] , s[x,y,0]))]



if __name__ == "__main__":
        f = open("average_misorientation_angle.txt", "w" )
        f.write("X\tY\ttheta \n")
        average_misorientation = np.zeros((r + 1, c + 1, 1))

        for x in range(1,r):
             for y in range(1,c,1):
                  for i in range(-1,2,1):
                        for j in range(-1,2,1):
                                if s[x,y,3] < 0.1 or s[x+i,y+j,3] < 0.1:
                                    average_misorientation[x,y,0] = average_misorientation[x,y,0] + 0
                                else:
                                     average_misorientation[x,y,0] = average_misorientation[x,y,0] + (theta(np.matmul(G[x,y,0],G[x+i,y+j,1])))
                                

                  average_misorientation[x,y,0] = (average_misorientation[x,y,0] - (theta(np.matmul(G[x,y,0],G[x,y,1]))))/8  
                                  
        for x in range(0,r+1):
             for y in range(0,c+1):
                  f.write("%s\t%s\t%s \n"%(x*stepsize,y*stepsize,np.degrees(average_misorientation[x,y,0])))