import numpy as np
import main
import random
import pandas as pd

n1 = 20  ### Number of X which should be used for avg 
n2 = 20 ### Number of Y which should be used for avg 


def average_grain_x(n1):
    random_y = random.sample(range(0,main.c),n1)
    subgrains =[]
    for y in random_y:
        i=0
        for x in range (0,main.r):
            #if main.s[x,y,3] < 0 or main.s[x+1,y,3] < 0:
                #pass
            if  np.degrees(main.theta(np.matmul(main.G[x,y,0],main.G[x+1,y,1]))) > 15:
                i+=1
        #print(y)        
        #print(i)          
        #print(main.r*main.stepsize_x/i)        
        subgrains.append(main.r*main.stepsize_x/i)
    
    #print(subgrains)
    #print(main.r*main.stepsize_x)    
    return sum(subgrains)/len(subgrains)
    
def average_grain_y(n2):
    random_x = random.sample(range(0,main.r),n2)
    subgrains =[]
    for x in random_x:
        i=0
        for y in range (0,main.c):
            #if main.s[x,y,3] < 0 or main.s[x,y+1,3] < 0:
            #   pass
            if np.degrees(main.theta(np.matmul(main.G[x,y,0],main.G[x,y+1,1]))) > 15:
                i+=1
        #print(x)        
        #print(i)        
        #print(main.r*main.stepsize_y/i) 
        subgrains.append(main.c*main.stepsize_y/i)
    #print(subgrains)    
    #print(main.r*main.stepsize_y)  
    return sum(subgrains)/len(subgrains)

if __name__ == "__main__":
    f = open("grain_size.txt", "w" )
    f.write("Avg Grain X = %s \n"%(average_grain_x(n1)))
    f.write("Avg Grain Y = %s \n"%(average_grain_y(n2)))