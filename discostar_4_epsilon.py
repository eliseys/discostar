import lmfit
import numpy as np
import discostar4_lib as d4l

h = 0.023
R = 0.23

phase = {
    '00_05.dat':[90,	 0,     0.0, 30.0,  0.0, 20.0],
#    '05_10.dat':[72,	 18,    0.0, 30.0,  0.0, 20.0],
#    '10_15.dat':[54,	 36,    0.0, 30.0,  0.0, 20.0],
#    '15_20.dat':[36,     54,    0.0, 30.0,  0.0, 20.0],
#    '20_25.dat':[18,     72,    0.0, 30.0,  0.0, 20.0],
#    '25_30.dat':[0,	 90,    0.0, 30.0,  0.0, 20.0],
#    '30_35.dat':[342,	 108,   0.0, 30.0,  0.0, 20.0],
#    '35_40.dat':[324,	 126,   0.0, 30.0,  0.0, 10.0],
#    '40_45.dat':[306,	 144,   0.0, 30.0,  -5.0, 5.0],
#    '45_50.dat':[288,	 162,   0.0, 30.0,  -10.0, 0.0],
#    '50_55.dat':[270,	 180,   0.0, 30.0,  -10.0, 0.0],
#    '55_60.dat':[252,	 198,   0.0, 30.0,  -10.0, 0.0],
#    '60_65.dat':[234,	 216,   0.0, 30.0,  -10.0, 0.0],
#    '65_70.dat':[216,	 234,   0.0, 30.0,  -10.0, 0.0],
#    '70_75.dat':[198,	 252,   0.0, 30.0,  -10.0, 0.0],
#    '75_80.dat':[180,	 270,   0.0, 30.0,  -5.0, 5.0],
#    '80_85.dat':[162,	 288,   0.0, 30.0,  -5.0, 5.0],
#    '85_90.dat':[144,	 306,   0.0, 30.0,   0.0, 20.0],
#    '90_95.dat':[126,	 324,   0.0, 30.0,   0.0, 20.0],
#    '95_00.dat':[108,	 342,   0.0, 30.0,  0.0, 20.0],
}


theta_list = [-3.0,]
kappa_list = [15.0,]


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
                x, data = np.loadtxt(data, dtype=('float'), usecols=(4, 1), unpack=True)
                
                pars.add_many(
                    ('border', 0.13, False),
                    ('Lx', 1.5e+37, True, 0.5e+37, 9.0e+37, None, 1.0e+37),
                    ('PSI_pr', PSI_pr, False, 0, 360, None),
                    ('z_tilt', z_tilt, False, 0.0, 360.0, None),
                    ('y_tilt', 15.0, True, 0.0, 40.0, None, 2.5),
                    #('y_tilt2', 15.0, True, 0.0, 30.0, None, 2.5),
                    #('epsilon_in', 0., True, epsilon_in_min, epsilon_in_max, None, 2.5),
                    #('epsilon_out', 0., True, epsilon_out_min, epsilon_out_max, None, 2.5),
                    ('z_tilt2', -54.0, True, -130.0, 0.0, None, 5.0),
                    ('T_disk', 15000.0, True, 10000.0, 55000.0, None),
                    ('delta', 0.5, True, 0.0, 1.0),
                    #('a', 6.35e11, False),
                    #('sigma', 5.6704e-5, False),                    
                    #('R', 0.23, False),
                    #('H', 0.023, False),

                )

                #pars.add('y_tilt', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                #pars.add('y_tilt2', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')

                pars.add('y_tilt2', expr='y_tilt*delta')

                #pars.add('alpha', expr='arccos(cos(y_tilt)*cos(y_tilt2)+sin(y_tilt)*sin(y_tilt2)*cos(z_tilt2))')

                
                
                #pars.add('T_disk', expr='((Lx/sigma)*(alpha/(8.0*pi*pi*R*R*a*a) + (H*a)/(14.0*pi*R*R*R*a*a*a)))**(1.0/4.0)')

                pars.add('disk_flux', expr='0')





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
                    #if ((X > border) and (X < border3)) or ((X > border4) and (X < 1.0 - border)):
                        x_short.append(x[i])
                        data_short.append(data[i])

                data_short = np.array(data_short)

                N = len(data_short)
        
                f = d4l.lc_disk(z_tilt, PSI_pr, Lx, y_tilt, y_tilt2, z_tilt2, disk_flux, T_disk, theta, kappa)

                #if np.sum(np.square(data_short - f(x_short)))/N < 0.07:
                print '{:8f}'.format(np.sum(np.square(data_short - f(x_short)))/N), '\t', '{:8f}'.format(z_tilt), '\t', '{:8f}'.format(y_tilt),  '\t','{:8f}'.format(y_tilt2), '\t', '{:8f}'.format(z_tilt2), '\t', '{:E}'.format(Lx), '\t', '{:8f}'.format(T_disk)
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
    
            print(lmfit.fit_report(result))
            print "kappa", "\t", kappa
            print "ns_theta", "\t", theta
            #ci = lmfit.conf_interval(mini, result)
            #lmfit.printfuncs.report_ci(ci)



