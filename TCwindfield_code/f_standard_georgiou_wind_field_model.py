# import math
import matplotlib.pyplot as plt
import numpy as np

class sgwfm:
    def f_standard_georgiou_wind_field_model(Rmax, dp, B, Vc, Lat):

        # Define basic inputs
        # Rmax = 30000  # unit: m
        # dp = 5000  # unit: Pa
        # B = 1.5
        # Vc = 5  # unit: m/s
        # Lat = 35  # Latitude in degree (°)

        # Corilois force
        # Lat = 35
        f = 2 * 7.292 * 10 ** (-5) * np.sin(Lat / 180 * np.pi)
        pc = 101000 - dp

        # Calculate wind field
        th = np.arange(np.pi / 180, 2 * np.pi + np.pi / 180, np.pi / 180)
        r = np.arange(2000, 1000000 + 1, 2000)  # unit: m
        r1 = r
        tt, rr = np.meshgrid(th, r)
        th = tt
        r = rr

        # Gourgiou wind field
        # vsc1=(B*dep./1.15*(rm./r).^B.*exp(-(rm./r).^B)+((c.*sin(pi/2-th)-r*f)/2).^2).^0.5+(c*sin(pi/2-th)-r*f)/2;  % Georgiou wind velocity;
        # vsc2=(B*dep./1.15*(rm./r).^B.*exp(-(rm./r).^B)+((r*f)/2).^2).^0.5+(-r*f)/2;  % Holland wind velocity;
        a1 = B * dp / 1.15 * (Rmax / r) ** B * np.exp(-(Rmax / r) ** B)
        vG_Gor = (a1 + ((Vc * np.sin(np.pi / 2 - th) - r * f) / 2) ** 2) ** 0.5 + (Vc * np.sin(np.pi / 2 - th) - r * f) / 2
        # Holland wind field
        vG_Hol = (a1 + ((r * f) / 2) ** 2) ** 0.5 + (-r * f) / 2


        plot_id = 0
        if plot_id == 1:
            # plot the wind field
            def polar_to_cartesian(r, theta):
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                return x, y


            xx, yy = polar_to_cartesian(r, th)

            # plt.figure(1, figsize=(8, 8))
            plt.figure(1)
            plt.contourf(xx / 1000, yy / 1000, vG_Gor)
            plt.xlabel("X (km)")
            plt.ylabel("Y (km)")
            plt.title('Giourgiou wind field')
            plt.colorbar(label='Wind speed (m/s)')
            plt.xlim(-4 * Rmax / 1000, 4 * Rmax / 1000)
            plt.ylim(-4 * Rmax / 1000, 4 * Rmax / 1000)


            plt.figure(2)
            # plt.figure(1, figsize=(6, 6))
            plt.contourf(xx / 1000, yy / 1000, vG_Hol)
            plt.xlabel("X (km)")
            plt.ylabel("Y (km)")
            plt.title('Holland wind field')
            plt.colorbar(label='Wind speed (m/s)')
            plt.xlim(-4 * Rmax / 1000, 4 * Rmax / 1000)
            plt.ylim(-4 * Rmax / 1000, 4 * Rmax / 1000)
            plt.show()

        return vG_Gor, vG_Hol, r1
