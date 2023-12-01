import esqg
import numpy as np
import matplotlib.pyplot as plt


# Create an instance of isqg_data
e = esqg_data()

# Define variables for the instance
N = 50
d.z = linspace(-500, -10, N)
d.c = 1 # given constant
d.n0 = 0.003  # given constant
d.lon = arange(0, 7, 0.1)
d.lat = arange(30, 37, 0.1)
d.rho0 = 1025.
d.window_width = 10
d.ssh = np.ones((d.lat.shape[0], d.lon.shape[0]))

# Solve the eSQG equation to obtain the density and velocity field
e = esqg.esqg_data()
e.solve_esqg()

# Plot the results
figure(0)
subplot(311)
plt.contourf(lon, lat, d.rho[-1, :, :])
plt.title('Surface density')
subplot(312)
plt.contourf(lon, lat, d.zeta[-1, :, :])
plt.title('Surface vorticity')
subplot(313)
plt.contourf(lon, lat, d.w[0, :, :])
plt.title('vertical velocity at 500 m')