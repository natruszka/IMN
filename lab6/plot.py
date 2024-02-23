import matplotlib.pyplot as plt
import numpy as np


file = "5a.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((50 + 1, 50 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 50, $\epsilon_1 = 1, \epsilon_2 = 1$")
plt.savefig("5a.png")
plt.clf()


file = "5b.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((100 + 1, 100 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 100, $\epsilon_1 = 1, \epsilon_2 = 1$")
plt.colorbar()
plt.savefig("5b.png")
plt.clf()


file = "5c.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((200 + 1, 200 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 200, $\epsilon_1 = 1, \epsilon_2 = 1$")
plt.savefig("5c.png")
plt.clf()


file = "6a.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((100 + 1, 100 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 100, $\epsilon_1 = 1, \epsilon_2 = 1$")
plt.savefig("6a.png")
plt.clf()


file = "6b.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((100 + 1, 100 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 100, $\epsilon_1 = 2, \epsilon_2 = 2$")
plt.colorbar()
plt.savefig("6b.png")
plt.clf()


file = "6c.txt"
data = np.loadtxt(file)
l = data[:,0]
V = data[:,1:]
V = V.reshape((100 + 1, 100 + 1))
x = np.arange(0, len(V)) * 0.1
y = np.arange(0, len(V[0])) * 0.1
plt.pcolor(x, y, V, cmap='bwr')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("nx = ny = 100, $\epsilon_1 = 10, \epsilon_2 = 10$")
plt.savefig("6c.png")
plt.clf()


print("done :)")
