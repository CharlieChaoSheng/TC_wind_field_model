
import numpy as np
# from f_standard_georgiou_wind_field_model import sgwfm
from f_georgiou_wind_field_model import gwfm
from scipy.interpolate import griddata

class ewugm:
    def estimate_wind_using_georgiou_model(lats, lons, lat, lon, Rmax, dp, B, Vc, Lat, heading, Factor_G2S, Factor_T):
        # # define basic inputs
        # lats = np.array([20, 20.4])
        # lons = np.array([120, 120.5])
        # lat = 19
        # lon = 119
        # Factor_G2S = 0.85
        # Factor_T = 0.93 # from the 1-h to 1 min
        #
        # Rmax = 30000  # unit: m
        # dp = 5000 # unit: Pa
        # B = 1.5
        # Vc = 5 # unit: m/s
        # Lat = 35 # Latitude in degree (°)
        # heading = 180 # degree 0-360 from true North

        # -----------------------------------------------------------------------------
        # calculate the angle and distance from the TC centre to the estimated sites
        a_earth = 6378.137 # km
        b_earth = 6356.752314 # km
        e2=(a_earth ** 2 - b_earth ** 2) / (a_earth ** 2)
        lat1 = lat * np.pi / 180
        dlat = np.pi * a_earth * (1-e2) / (180 * (1-e2 * (np.sin(lat1))**2)**(3/2))
        dlon = np.pi * a_earth * np.cos(lat1) / (180 * (1-e2 * (np.sin(lat1))**2)**(1/2))

        pdistx = (lons - lon) * dlon * 1000 # m
        pdisty = (lats - lat) * dlat * 1000 # m
        pdist = np.sqrt(pdistx**2 + pdisty**2)

        ind = np.where(pdist <= 250000) # signal the sites within 250km of typhoon center
        pdistx = pdistx[ind]
        pdisty = pdisty[ind]

        bm = np.zeros(pdisty.shape)
        bm[np.where(pdisty < 0)] = 1
        pangle = (-1)**bm * np.arccos(pdistx / pdist) * 180 / np.pi
        pangle[np.logical_and(pangle > 179, pangle <= 180)] = 179
        pangle[np.isnan(pangle)] = 0

        # obtain gorgiou wind field
        G_Gor, thetat1, r1 = gwfm.f_georgiou_wind_field_model(Rmax, dp, B, Vc, Lat, heading)
        r2, thetat2 = np.meshgrid(r1, thetat1)

        points = np.column_stack((r2.flatten(), thetat2.flatten()))
        values = G_Gor.flatten()

        # VG1 = griddata(points, values, points, method="linear")
        # V_check = VG1 - values
        VG1 = griddata(points, values, (pdist, pangle), method="linear")

        VG1[np.isnan(VG1)] = 0
        VG1[pdist == 0] = 0
        # VG1[np.where(pdist == 0)] = 0
        V_Gor_G = VG1

        # Conversion From Gradient height to surface
        V_Gor_S = V_Gor_G * Factor_G2S
        # Conversion for different time duration
        V_Gor_S_T = V_Gor_S * Factor_T

        # output V_Gor_S_T
        return V_Gor_S_T

