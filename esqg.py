class esqg_data:
    """
    Class for reconstructing a 3D ocean field using surface observations.
    Initialization requires effective buoyancy frequency N0, longitude and latitude, 
    surface height, and the amplitude of the vertical velocity.

    Args:
        n0 (float): effective buoyancy frequency
        c (float): amplitude of vertical velocity
        z (float): z value
        ssh (float): sea surface height
        lon (float): longitude
        lat (float): latitude
        rho0 (float): mean density
        window_width (int): the width of trapezoid window

    """

    def __init__(self, n0: float = 0.0, c: float = 0.0, z: float = 0.0, ssh: float = 0.0,
                 lon: float = 0.0, lat: float = 0.0, rho0: float = 1025.0, window_width: int = 10) -> None:
        """
        Args:
            n0 (float): effective buoyancy frequency
            c (float): amplitude of vertical velocity
            z (float): z value
            ssh (float): sea surface height
            lon (float): longitude
            lat (float): latitude
            rho0 (float): mean density
            window_width (int): the width of trapezoid window
        Returns:
            None
        """

        # Set attributes
        self.n0 = n0
        self.c = c
        self.z = z
        self.rho0 = rho0
        self.ssh = ssh
        self.lon = lon
        self.lat = lat
        self.window_width = window_width
        return

    def precondition(self):
        """
        Set up vertical difference matrix, initialize parameters and check N2.

        :return: None
        """
        import numpy as np
        from numpy import fft
        import seawater as sw

        # Set coordinate
        self.surface = np.argsort(abs(self.z))[0]
        if np.isscalar(self.lon) and np.isscalar(self.x):
            print("error:  please specify the latitude and longitude")
        elif not np.isscalar(self.lon):
            self.x, self.y = lonlat2xy(self.lon, self.lat)
        else:
            self.x, self.y = np.meshgrid(self.x, self.y)

        self.dx = (self.x[:, 1] - self.x[:, 0]).reshape(-1, 1)
        self.dy = (self.y[1, :] - self.y[0, :]).reshape(1, -1)
        self.ny, self.nx = self.x.shape

        # Initialize parameters
        self.f0 = lat2f(self.lat.mean())
        self.g = 9.81

        # Set wavenumbers
        ny, nx = self.ssh.shape
        k = np.arange(-nx / 2, nx / 2, 1)
        kx = 2 * np.pi * k[None, :] / nx / self.dx
        k = np.arange(-ny / 2, ny / 2, 1)
        ky = 2 * np.pi * k[:, None] / ny / self.dy
        self.wv2 = (kx ** 2 + ky ** 2) ** 0.5

        # Apply 2D trapezoid window
        num = self.window_width
        winx = np.ones((self.nx))
        winx[0:num] = np.linspace(0, 1, num)
        winx[-num::] = np.linspace(1, 0, num)
        winy = np.ones((self.ny))
        winy[0:num] = np.linspace(0, 1, num)
        winy[-num::] = np.linspace(1, 0, num)
        wx, wy = np.meshgrid(winx, winy)
        self.window = wx * wy
        self.ssh = self.ssh * self.window
        
        return

    def solve_esqg(self):
        import numpy as np
        from numpy import fft

        self.precondition()

        N0 = self.n0
        k = self.wv2
        z = self.z
        c = self.c
        f0 = self.f0
        self.ddy, self.ddx = np.meshgrid(self.dy[0, :], self.dx[:, 0])
        ##################################################
        sshhat = fft.fftshift(fft.fft2(self.ssh))
        psihat = self.g / f0 * sshhat[None, :, :] * np.exp(N0 / f0 * k[None, :, :] * z[:, None, None])
        bhat = N0 * k[None, :, :] / c * psihat
        self.psi = fft.ifft2(fft.ifftshift(psihat, axes = (1, 2)), axes = (1, 2)).real
        self.b = fft.ifft2(fft.ifftshift(bhat, axes = (1, 2)), axes = (1, 2)).real
        
        self.rho = -self.b / self.g * self.rho0
        self.u, self.v = psi2uv(self.lon, self.lat, self.psi)

        vorthat = - (k[None, :, :] ** 2) * psihat
        self.vort = fft.ifft2(fft.ifftshift(vorthat, axes = (1, 2)), axes = (1, 2)).real
        ##################################################
        psix = (self.psi - np.roll(self.psi, 1, axis = -1)) / self.ddx[None, :, :]
        psiy = (self.psi - np.roll(self.psi, 1, axis = 1)) / self.ddy[None, :, :]
        bx = (self.b - np.roll(self.b, 1, axis = -1)) / self.ddx[None, :, :]
        by = (self.b - np.roll(self.b, 1, axis = 1)) / self.ddy[None, :, :]
        jz = fft.fftshift(fft.fft2(psix * by - psiy * bx, axes = (1, 2)), axes = (1, 2))
        js = jz[self.surface, :, :]
        what = - c ** 2 / (N0 ** 2) * (-js[None, :, :] * np.exp(N0 / f0 * k[None, :, :] * z[:, None, None]) + jz)
        self.w = fft.ifft2(fft.ifftshift(what, axes = (1, 2)), axes = (1, 2)).real
        return

def lonlat2xy(lon,lat):
    """
    Converts the given longitude and latitude arrays to x and y coordinates
    using the WGS84 reference ellipsoid.

    Parameters:
        lon (ndarray): Array of longitudes
        lat (ndarray): Array of latitudes

    Returns:
        x (ndarray): Array of x-coordinates
        y (ndarray): Array of y-coordinates
    """
    import numpy as np
    r = 6371.e3
    lon = lon-lon[0]
    if lon.ndim == 1:
        lon,lat = np.meshgrid(lon,lat)
    x = 2*np.pi*r*np.cos(lat*np.pi/180.)*lon/360.
    y = 2*np.pi*r*lat/360.
    return x,y

def lat2f(d):
    """
    Calculates the Coriolis parameter from the latitude d.

    Parameters:
        d (ndarray): Array of latitudes (1- or 2-D)

    Returns:
        f (ndarray): Array of Coriolis parameter values
    """
    from math import pi, sin
    return 2.0*0.729e-4*sin(d*pi/180.0)

def gradxy(lon,lat,d):
    """
    Calculate the x and y gradients of a field
    lon: longitudes, 1- or 2-D
    lat: latitudes, 1- or 2-D
    d: field, 2- or 3-D

    The calculated velocity fields share the same dimensions with d,
    i.e. the velocity fields are shifted.
    """
    from numpy import c_,r_,diff,zeros
    x,y = lonlat2xy(lon,lat)
    def uv(x,y,d):
        x = c_[x, x[:,-2]]
        y = r_[y, y[-2,:].reshape(1,-1)]

        sshx = c_[d, d[:,-1]]
        v = diff(sshx,axis=1)/diff(x,axis=1)
        v[:,-1] = v[:,-2]

        sshy = r_[d, d[-1,:].reshape(1,-1)]
        u = diff(sshy,axis=0)/diff(y,axis=0)
        u[-2,:]=u[-1,:]
        return u,v
    if d.ndim == 2:
        return uv(x,y,d)
    elif d.ndim == 3:
        u, v = zeros(d.shape), zeros(d.shape)
        for i in range(d.shape[0]):
            ut,vt = uv(x,y,d[i,:,:])
            v[i,:,:] = vt
            u[i,:,:] = ut
        return u, v

def psi2uv(lon,lat,d):
    """
    Calculate horizontal velocity from streamfunction.

    Parameters
    ----------
    lon : 1D or 2D array-like
        Longitudes
    lat : 1D or 2D array-like
        Latitudes
    d : 2D or 3D array-like
        Streamfunction

    Returns
    -------
    u : 2D or 3D array-like
        Zonal velocity
    v : 2D or 3D array-like
        Meridional velocity
    """
    u, v = gradxy(lon, lat, d)
    u = -1 * u
    return u, v


def remove_largescale(data):
    """
    remove the large-scale background signals through the bilinear least-squares fitting

    Parameters
    ----------
    data : 2D or 3D array-like
        background information

    Returns
    -------
    data : same as input but removing the large-scale signals
    data_mean: the large-scale signals
    """
    from sklearn import linear_model
    import numpy as np
    data[np.isnan(data)] = np.nanmean(data)
    if data.ndim == 2:
        xx, yy = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))
        X, Z = np.column_stack((xx.flatten(), yy.flatten())), data.flatten()
        regr = linear_model.LinearRegression()
        regr.fit(X, Z)
        a, b = regr.coef_, regr.intercept_
        data_mean = regr.predict(X).reshape(np.shape(data))
        data = data - data_mean
    else:
        x = []
        for i in range(data.shape[0]):
            xx, yy = np.meshgrid(np.arange(data[i].shape[0]), np.arange(data[i].shape[1]))
            X, Z = np.column_stack((xx.flatten(), yy.flatten())), data[i].flatten()
            regr = linear_model.LinearRegression()
            regr.fit(X, Z)
            a, b = regr.coef_, regr.intercept_
            data_mean = regr.predict(X).reshape(np.shape(data[i]))
            x.append(data_mean)
        data_mean = np.array(x)
        data = data - data_mean
    return data, data_mean
