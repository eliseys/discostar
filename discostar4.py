import lmfit
import numpy as np
import discostar4_lib as d4l

phase = {
#    '25_30.dat':[0,	270],
#    '20_25.dat':[342,	288],
#    '15_20.dat':[324,	306],
    '10_15.dat':[306,	324],
#    '05_10.dat':[288,	342],
#    '00_05.dat':[270,	0],
#    '95_00.dat':[252,	18],
#    '90_95.dat':[234,	36],
#    '85_90.dat':[216,	54],
#    '80_85.dat':[198,	72],
#    '75_80.dat':[180,	90],
#    '70_75.dat':[162,	108],
#    '65_70.dat':[144,	126],
#    '60_65.dat':[126,	144],
#    '55_60.dat':[108,	162],
#    '50_55.dat':[90,	180],
#    '45_50.dat':[72,    198],
#    '40_45.dat':[54,	216],
#    '35_40.dat':[36,	234],
#    '30_35.dat':[18,	252]
}

theta_list = [0,]
kappa_list = [-15,]

for theta in theta_list:
    for kappa in kappa_list:
        for i, x in enumerate(phase):

            z_tilt = phase[x][0]
            PSI_pr = phase[x][1]
            data_file = x
            
            print  z_tilt, PSI_pr, data_file
            
            path_to_data_file = './OBSERVED_DATA/'+data_file
            
            pars = lmfit.Parameters()
            
            with open(path_to_data_file, 'r') as data:
                x, data = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)
                
                pars.add_many(
                    ('border', 0.25, False),
                    ('Lx', 2.0e+37, True, 1.0e+37, 9.0e+37, None, 1.0e+37),
                    ('PSI_pr', PSI_pr, False),
                    ('z_tilt', z_tilt, False),
                    #('y_tilt', 4.27259339, True, 0.0, 25.0, None, 2.5),
                    ('y_tilt2', 4.27259339, True, 0.0, 25.0, None, 2.5),
                    ('z_tilt2', 0.0, True, -90.0, 90.0, None, 5.0),
                    ('T_disk', 28175.4064, True, 10000.0, 80000.0, None)
                )

                pars.add('y_tilt', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                
                #pars.add('y_tilt2', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                pars.add('disk_flux', expr='0')

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
        
                f = d4l.lc_disk(z_tilt, PSI_pr, Lx, y_tilt, y_tilt2, z_tilt2, disk_flux, T_disk, theta, kappa)

                #if np.sum(np.square(data_short - f(x_short)))/N < 0.07:
                print np.sum(np.square(data_short - f(x_short)))/N, y_tilt, y_tilt2, z_tilt2, T_disk, Lx, disk_flux
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
            print "kappa", "\t", kappa
            print "ns_theta", "\t", theta
            #ci = lmfit.conf_interval(mini, result)
            #lmfit.printfuncs.report_ci(ci)







#def residual(pars, x, data):
#
#    z_tilt = pars['z_tilt'].value
#    PSI_pr = pars['PSI_pr'].value
    
#    Lx = pars['Lx'].value
#    y_tilt2 = pars['y_tilt2'].value
#    z_tilt2 = pars['z_tilt2'].value

#    y_tilt = pars['y_tilt'].value
#    T_disk = pars['T_disk'].value
#    
#    f = d4l.lc(Lx, y_tilt2, z_tilt2, y_tilt, T_disk, z_tilt, PSI_pr)

#    return data - f(x)


#min_ = lmfit.Minimizer(residual, pars, fcn_args=(x, data))
#result = min_.minimize()



#print '\n'


#print result.params.pretty_print

#print(lmfit.fit_report(result.params))

#ci = lmfit.conf_interval(min_, result)
#lmfit.printfuncs.report_ci(ci)
