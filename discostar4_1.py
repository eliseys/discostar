import lmfit
import numpy as np
import discostar4_lib as d4l

phase = {
    1:['75_80.dat', 180, 90],
    2:['80_85.dat', 198, 72],
    3:['85_90.dat', 216, 54],
    4:['90_95.dat', 234, 36],
    5:['95_00.dat', 252, 18],
    6:['00_05.dat', 270, 0],
    7:['05_10.dat', 288, 342],
    8:['10_15.dat', 306, 324],
    9:['15_20.dat', 324, 306],
    10:['20_25.dat', 342, 288],
    11:['25_30.dat', 0, 270],
    12:['30_35.dat', 18, 252],
    13:['35_40.dat', 36, 234],
    14:['40_45.dat', 54, 216],
    15:['45_50.dat', 72, 198],
    16:['50_55.dat', 90, 180],
    17:['55_60.dat', 108, 162],
    18:['60_65.dat', 126, 144],
    19:['65_70.dat', 144, 126],
    20:['70_75.dat', 162, 108]
}

Lx_in = 3.6029e+37
y_tilt_in = 19.7117254
y_tilt2_in = 8.12772819
z_tilt2_in = 89.7942076
T_disk_in = 13332.7143

delta = 0.3

for i in range(1,21):

    data_file = phase[i][0]
    z_tilt = phase[i][1]
    PSI_pr = phase[i][2]

    print  z_tilt, PSI_pr, data_file

    path_to_data_file = './OBSERVED_DATA/'+data_file

    pars = lmfit.Parameters()
    
    with open(path_to_data_file, 'r') as data:
        x, data = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)

        pars.add_many(
            ('border', 0.13, False),
            ('Lx', 6.23e+37, True, Lx_in*(1.0 - delta), Lx_in*(1.0 + delta), None, 1.0e+37),
            ('PSI_pr', PSI_pr, False),
            ('z_tilt', z_tilt, False),
            ('y_tilt', 14.0, True, y_tilt_in*(1.0 - delta), y_tilt_in*(1.0 + delta), None, 2.5),
            ('y_tilt2', 14.7, True, y_tilt2_in*(1.0 - delta), y_tilt2_in*(1.0 + delta), None, 2.5),
            ('z_tilt2', 52.0, True, z_tilt2_in*(1.0 - delta), z_tilt2_in*(1.0 + delta), None, 5.0),
            ('T_disk', 15000.0, True, T_disk_in*(1.0 - delta), T_disk_in*(1.0 + delta), None)
        )

        pars.add('disk_flux', expr='(5.79414e+38 * abs(cos(y_tilt*pi/180.0 - 80.0*pi/180.0) * cos(z_tilt)) + 8.50668e+37) * (T_disk/10000.0)**4.0')

        #('disk_flux', 38.079181246047625, True, 25.0, 39.477121254719663, None, 1.4477121254719663)
        #a               = 5.79414e+38      +/- 1.475e+36    (0.2545%)
        #b               = 8.50668e+37      +/- 6.896e+35    (0.8106%)
        #('disk_flux', 1.2e+38, True, 1.0e+25, 3e+39, expr='5.79414e+38*(cos(y_tilt*pi/180.0 - 80.0*pi/180.0)*cos(z_tilt)) + 8.50668e+37')




        
    def test_fit(pars, x, data):
        print pars['y_tilt']
        return (pars['y_tilt']-5)**4 + 10000

    def residual(pars, x, data):

        z_tilt = pars['z_tilt'].value
        PSI_pr = pars['PSI_pr'].value

        Lx = pars['Lx'].value
        
        y_tilt = pars['y_tilt'].value
        y_tilt2 = pars['y_tilt2'].value
        z_tilt2 = pars['z_tilt2'].value

        T_disk = pars['T_disk'].value
        disk_flux = pars['disk_flux'].value

        border = pars['border'].value

        x_short = []
        data_short = []
    
        for i, X in enumerate(x):
            if (X > border) or (X < 1.0 - border):
                x_short.append(x[i])
                data_short.append(data[i])

        data_short = np.array(data_short)

        N = len(data_short)
        
        f = d4l.lc_disk(z_tilt, PSI_pr, Lx, y_tilt, y_tilt2, z_tilt2, disk_flux)

        #if np.sum(np.square(data_short - f(x_short)))/N < 0.07:
        #print np.sum(np.square(data_short - f(x_short)))/N, y_tilt, y_tilt2, z_tilt2, T_disk, Lx, disk_flux
        #print np.sum(np.square(data_short - f(x_short)))/N
        
        return data_short - f(x_short)



    mini = lmfit.Minimizer(residual, pars, fcn_args=(x, data))
    #result = mini.minimize(epsfcn=1.0e-4)
    result = mini.minimize(method='differential_evolution')
    #result = mini.minimize(method='brute')


    #mini = lmfit.Minimizer(test_fit, pars, fcn_args=(x, data))
    #result = mini.minimize()


    #result = mini.minimize(method='differential_evolution')
    
    #print result.params

    #print result.params.pretty_print

    print(lmfit.fit_report(result.params))

    #ci = lmfit.conf_interval(mini, result)
    #lmfit.printfuncs.report_ci(ci)

    #print result.params['Lx'].value, result.params['T_disk'].value
    
    Lx_in = result.params['Lx'].value
    y_tilt_in = result.params['y_tilt'].value
    y_tilt2_in = result.params['y_tilt2'].value
    z_tilt2_in = result.params['z_tilt2'].value
    T_disk_in = result.params['T_disk'].value
