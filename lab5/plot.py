import matplotlib.pyplot as plt
import numpy as np

delta = 0.2
nx = 128
ny = 128
xmax = delta * nx
ymax = delta * ny

file = "k16.txt"
data = np.loadtxt(file)

data = np.transpose(data)
x = np.arange(0, xmax + 0.1, 16 * delta)
plt.pcolor(x, x, data, cmap="jet") #nie obliczam y bo sa takie same
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("k=16")
plt.savefig("k16.png")
plt.clf()

file = "k8.txt"
data = np.loadtxt(file)

x = np.arange(0, xmax + 0.1, 8 * delta)
data = data.reshape(len(x), len(x))
data = np.transpose(data)
plt.pcolor(x, x, data, cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("k=8")
plt.savefig("k8.png")
plt.clf()

file = "k4.txt"
data = np.loadtxt(file)

data = np.transpose(data)
x = np.arange(0, xmax + 0.1, 4 * delta)
plt.pcolor(x, x, data, cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("k=4")
plt.savefig("k4.png")
plt.clf()

file = "k2.txt"
data = np.loadtxt(file)

data = np.transpose(data)
x = np.arange(0, xmax + 0.1, 2 * delta)
plt.pcolor(x, x, data, cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("k=2")
plt.savefig("k2.png")
plt.clf()

file = "k1.txt"
data = np.loadtxt(file)

data = np.transpose(data)
x = np.arange(0, xmax + 0.1, 1 * delta)
plt.pcolor(x, x, data, cmap="jet")
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("k=1")
plt.savefig("k1.png")
plt.clf()

file = "it0.txt"
data = np.loadtxt(file)
y = data[:, 0]
x = data[:,1]
plt.plot(x,y, label = "k=16")

file = "it1.txt"
data = np.loadtxt(file)
y = data[:, 0]
x = data[:,1]
plt.plot(x,y, label = "k=8")

file = "it2.txt"
data = np.loadtxt(file)
y = data[:, 0]
x = data[:,1]
plt.plot(x,y, label = "k=4")

file = "it3.txt"
data = np.loadtxt(file)
y = data[:, 0]
x = data[:,1]
plt.plot(x,y, label = "k=2")

file = "it4.txt"
data = np.loadtxt(file)
y = data[:, 0]
x = data[:,1]
plt.plot(x,y, label = "k=1")
plt.legend()
plt.title("S(it)")
plt.grid()
plt.savefig("it.png")

