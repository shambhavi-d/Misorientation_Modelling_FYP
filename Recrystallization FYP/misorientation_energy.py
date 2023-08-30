import pandas as pd
import os
import numpy as np

path = os.getcwd() + "\Recrystallization FYP\HR.xlsx"  #D:\Python Codes\Recrystallization FYP\HR.xlsx
#print(path)
df = pd.read_excel(path)
# print(df)
df = df.to_numpy()

stepsize = df[1,3] - df[0,3]

df[:,3] = (df[:, 3] / (stepsize)).astype(int)
df[:,4] = (df[:, 4] / (stepsize)).astype(int)    ##### makes the step size to be one (for x,y co-ordinates doesnt affect euler angles) Needs to be generalized

r,c = int(np.max(df[:,3])),int(np.max(df[:,4]))
# print(r,c)
s = np.zeros((r+1,c+1,3))
for i in df:
    # print(i)
    s[int(i[3])][int(i[4])] = i[0:3]

def g(phi1,phi,phi2):
    g_one = np.array([ [np.cos(np.degrees(phi1)), np.sin(np.degrees(phi1)), 0] ,[-np.sin(np.degrees(phi1)), np.cos(np.degrees(phi1)), 0] ,[0,0,1]])
    g_two = np.array([ [np.cos(np.degrees(phi2)), np.sin(np.degrees(phi2)), 0] ,[-np.sin(np.degrees(phi2)), np.cos(np.degrees(phi2)), 0] ,[0,0,1]])
    g = np.array([ [1,0,0], [0, np.cos(np.degrees(phi)), np.sin(np.degrees(phi))] ,[0, -np.sin(np.degrees(phi)), np.cos(np.degrees(phi))]])

    b = np.matmul(g,g_one)
    c = np.matmul(g_two,b)

    return c

G = np.zeros((r+1,c+1,2))



print(g(s[1,1,0],s[1,1,1],s[1,1,2]))