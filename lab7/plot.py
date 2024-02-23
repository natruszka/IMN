import matplotlib.pyplot as plt
import numpy as np

file = "psi-1000.txt"
psi = np.loadtxt(file).reshape((201, 91))
psi = np.where(psi==0, np.nan, psi)
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi),  levels = 50)
plt.colorbar()
plt.title("psi -1000")
plt.savefig("psi-1000.png")
plt.clf()

file = "zeta-1000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi),  levels = 50)
plt.colorbar()
plt.title("zeta -1000")
plt.savefig("zeta-1000.png")
plt.clf()

file = "v-1000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("v -1000")
plt.savefig("v-1000.png")
plt.clf()

file = "u-1000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("u -1000")
plt.savefig("u-1000.png")
plt.clf()

file = "psi4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
psi = np.where(psi==0, np.nan, psi)
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi),  levels = 50)
plt.colorbar()
plt.title("psi 4000")
plt.savefig("psi4000.png")
plt.clf()

file = "zeta4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi), levels = 50)
plt.colorbar()
plt.title("zeta 4000")
plt.savefig("zeta4000.png")
plt.clf()

file = "v4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("v 4000")
plt.savefig("v4000.png")
plt.clf()

file = "u4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("u 4000")
plt.savefig("u4000.png")
plt.clf()

file = "psi-4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
psi = np.where(psi==0, np.nan, psi)
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi), levels = 50)
plt.colorbar()
plt.title("psi -4000")
plt.savefig("psi-4000.png")
plt.clf()

file = "zeta-4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.contour(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi),  levels = 50)
plt.colorbar()
plt.title("zeta -4000")
plt.savefig("zeta-4000.png")
plt.clf()

file = "v-4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("v -4000")
plt.savefig("v-4000.png")
plt.clf()

file = "u-4000.txt"
psi = np.loadtxt(file).reshape((201, 91))
x = np.arange(0, 200 * 0.01 + 0.01, 0.01)
y = np.arange(0, 90 * 0.01 + 0.01, 0.01)
plt.pcolor(x, y, np.transpose(psi), vmin=np.nanmin(psi), vmax=np.nanmax(psi))
plt.colorbar()
plt.title("u -4000")
plt.savefig("u-4000.png")
plt.clf()

