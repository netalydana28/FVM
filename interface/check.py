import numpy as np
x = [170, 175, 160, 165]

for i in range(len(x)):
    print(f"degrees: {x[i]} cos: {np.cos(x[i]* np.pi/180)} sin: {np.sin(x[i]* np.pi/180)}\n")