
# generate the standard georgiou wind field and rotate with translation direction

import numpy as np
from f_standard_georgiou_wind_field_model import sgwfm

class gwfm:
    def f_georgiou_wind_field_model(Rmax, dp, B, Vc, Lat, heading):
        # Rmax = 30000  # unit: m
        # dp = 5000 # unit: Pa
        # B = 1.5
        # Vc = 5 # unit: m/s
        # Lat = 35 # Latitude in degree (°)
        # heading = 60 # degree 0-360 from true North

        G_Gor, G_Hol, r1 = sgwfm.f_standard_georgiou_wind_field_model(Rmax, dp, B, Vc, Lat)
        G_Gor = G_Gor.T
        G_Hol = G_Hol.T

        theta1 = np.arange(1, 360+0.001).T - heading
        theta1_O = theta1.T

        theta1 = np.where((theta1 > 0) & (theta1 < 360), np.mod(theta1, 180) + (theta1 > 180) * (-180), theta1)
        theta1 = np.where(theta1 > 360, np.mod(theta1, 360), theta1)
        theta1 = np.where((theta1 < 0) & (theta1 > -360), np.mod(theta1, -180) + (theta1 < -180) * 180, theta1)
        theta1 = np.where(theta1 < -360, np.mod(theta1, -360), theta1)
        theta1 = np.where(theta1 == -180, -180, theta1)
        theta1 = np.where(theta1 == 180, -180, theta1)

        tem1 = np.column_stack((theta1, G_Gor))
        tem2 = np.column_stack((theta1, G_Hol))
        t1 = tem1[:, 0].argsort()
        tem1 = tem1[tem1[:, 0].argsort()]
        tem2 = tem2[tem2[:, 0].argsort()]
        thetat1 = tem1[:, 0]

        G_Gor = tem1[:, 1:]
        G_Hol = tem2[:, 1:]

        plot_id = 1
        if plot_id == 1:
            import matplotlib.pyplot as plt

            thetat1_plot = np.append(thetat1, [thetat1[0] + 360])
            th = thetat1_plot * np.pi / 180
            [tt, rr] = np.meshgrid(th, r1)

            # plot the wind field
            def polar_to_cartesian(r, theta):
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                return x, y

            xx, yy = polar_to_cartesian(rr, tt)

            G_Gor_plot = np.row_stack((G_Gor, G_Gor[0, :]))
            G_Hol_plot = np.row_stack((G_Hol, G_Hol[0, :]))
            G_Gor_plot = G_Gor_plot.T
            G_Hol_plot = G_Hol_plot.T

            # import pandas as pd
            # data_df = pd.DataFrame(G_Gor_plot.T)
            # writer = pd.ExcelWriter('Check_WF.xlsx')
            # data_df.to_excel(writer,'page_1', float_format = '%0.5f')
            # writer.save()


            # plt.figure(1, figsize=(8, 8))
            fig = plt.figure(1)
            cp = plt.contourf(xx / Rmax, yy / Rmax, G_Gor_plot)
            # levels = np.arange(0, 55 + 1, 5)
            # plt.clabel(cp, inline=True, fontsize=10)
            plt.xlabel("x/Rmax")
            plt.ylabel("y/Rmax")
            plt.title('Giourgiou wind field')
            plt.colorbar(label='Wind speed (m/s)')
            plt.xlim(-4, 4)
            plt.ylim(-4, 4)
            # plt.show()


            plt.figure(2)
            # plt.figure(1, figsize=(6, 6))
            plt.contourf(xx / 1000, yy / 1000, G_Hol_plot)
            plt.xlabel("x/Rmax ")
            plt.ylabel("y/Rmax")
            plt.title('Holland wind field')
            plt.colorbar(label='Wind speed (m/s)')
            plt.xlim(-4 * Rmax / 1000, 4 * Rmax / 1000)
            plt.ylim(-4 * Rmax / 1000, 4 * Rmax / 1000)
            plt.show()

        cs = 1
        return G_Gor, thetat1, r1



