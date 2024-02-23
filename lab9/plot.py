import matplotlib.pyplot as plt
import numpy as np
values = [100,200,500,1000,2000]
for value in values:
    file = "t" + str(value) +".txt"
    t = np.loadtxt(file).reshape((41, 41))
    x = np.arange(0, 41)
    y = np.arange(0, 41)
    plt.pcolor(x, y, np.transpose(t), cmap="inferno")
    plt.colorbar()
    plt.title("it = " + str(value))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("t"+str(value)+".png")
    plt.clf()
    file = "dt" + str(value) +".txt"
    t = np.loadtxt(file).reshape((41, 41))
    x = np.arange(0, 41)
    y = np.arange(0, 41)
    plt.pcolor(x, y, np.transpose(t), cmap="inferno")
    plt.colorbar()
    plt.title("it = " + str(value))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("dt"+str(value)+".png")
    plt.clf()