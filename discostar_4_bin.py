import lmfit
import numpy as np
import discostar4_lib as d4l

h = 0.023
R = 0.23

#################################################################
# BIN BIN BIN
#################################################################


phase = {
    'bin20__00_05.dat':[90,	 0,     0.0, 30.0,  0.0, 20.0],
    'bin20__05_10.dat':[72,	 18,    0.0, 30.0,  0.0, 20.0],
    'bin20__10_15.dat':[54,	 36,    0.0, 30.0,  0.0, 20.0],
    'bin20__15_20.dat':[36,      54,    0.0, 30.0,  0.0, 20.0],
    'bin20__20_25.dat':[18,	 72,    0.0, 30.0,  0.0, 20.0],
    'bin20__25_30.dat':[0,	 90,    0.0, 30.0,  0.0, 20.0],
    'bin20__30_35.dat':[342,	 108,   0.0, 30.0,  0.0, 20.0],
    'bin20__35_40.dat':[324,	 126,   0.0, 30.0,  0.0, 10.0],
    'bin20__40_45.dat':[306,	 144,   0.0, 30.0,  -5.0, 5.0],
    'bin20__45_50.dat':[288,	 162,   0.0, 30.0,  -10.0, 0.0],
    'bin20__50_55.dat':[270,	 180,   0.0, 30.0,  -10.0, 0.0],
    'bin20__55_60.dat':[252,	 198,   0.0, 30.0,  -10.0, 0.0],
    'bin20__60_65.dat':[234,	 216,   0.0, 30.0,  -10.0, 0.0],
    'bin20__65_70.dat':[216,	 234,   0.0, 30.0,  -10.0, 0.0],
    'bin20__70_75.dat':[198,	 252,   0.0, 30.0,  -10.0, 0.0],
    'bin20__75_80.dat':[180,	 270,   0.0, 30.0,  -5.0, 5.0],
    'bin20__80_85.dat':[162,	 288,   0.0, 30.0,  -5.0, 5.0],
    'bin20__85_90.dat':[144,	 306,   0.0, 30.0,   0.0, 20.0],
    'bin20__90_95.dat':[126,	 324,   0.0, 30.0,   0.0, 20.0],
    'bin20__95_00.dat':[108,	 342,   0.0, 30.0,  0.0, 20.0],
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
                x, data, weight = np.loadtxt(data, dtype=('float'), usecols=(0,1,2), unpack=True)
                
                pars.add_many(
                    ('border', 0.13, False),
                    ('Lx', 2.0e+37, True, 0.5e+37, 5.0e+37, None, 1.0e+37),
                    ('PSI_pr', PSI_pr, False, 0, 360, None),
                    ('z_tilt', z_tilt, False, 0.0, 360.0, None),
                    ('y_tilt', 15.0, True, 0.0, 30.0, None, 2.5),
                    #('y_tilt2', 15.0, True, 0.0, 30.0, None, 2.5),
                    #('epsilon_in', 0., True, epsilon_in_min, epsilon_in_max, None, 2.5),
                    #('epsilon_out', 0., True, epsilon_out_min, epsilon_out_max, None, 2.5),
                    ('z_tilt2', -54.0, False, -110.0, 0.0, None, 5.0),
                    ('T_disk', 15000.0, True, 10000.0, 40000.0, None),
                    ('delta', 0.5, True, 0.0, 1.0)
                    #('inc', inc, False, 0, 90),
                )

                #pars.add('y_tilt', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                #pars.add('y_tilt2', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                pars.add('y_tilt2', expr='y_tilt*delta')


                #pars.add('cos_theta_minus_out', expr='(cos(inc)*sin(epsilon_out*(pi/180.0)) - sin(inc)*cos(z_tilt*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_out*(pi/180.0))**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0)')
                #pars.add('cos_theta_plus_out', expr='(cos(inc)*sin(epsilon_out*(pi/180.0)) + sin(inc)*cos(z_tilt*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_out*(pi/180.0))**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0)')
                #pars.add('cos_theta_minus_in', expr='(cos(inc)*sin(epsilon_in*(pi/180.0)) - sin(inc)*cos((z_tilt - z_tilt2)*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_in*(pi/180.0))**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0)')
                #pars.add('cos_theta_plus_in', expr='(cos(inc)*sin(epsilon_in*(pi/180.0)) + sin(inc)*cos((z_tilt - z_tilt2)*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_in*(pi/180.0))**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0)')

                #pars.add('y_tilt', expr='min((180.0/pi)*acos(cos_theta_minus_out), (180.0/pi)*acos(cos_theta_plus_out))')

                #pars.add('y_tilt2', expr='min((180.0/pi)*acos(cos_theta_minus_in), (180.0/pi)*acos(cos_theta_plus_in))')

                
                #pars.add('y_tilt', expr='-(180.0/pi)*asin((-cos(i)*sin(epsilon_out) + sin(i)*cos(z_tilt*(pi/180.0))*sqrt((cos(i)**2.0 - sin(epsilon_out)**2.0) + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0))/(cos(i)**2.0 + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0))')
                #pars.add('y_tilt', expr='sqrt((cos(i)**2.0 - sin(epsilon_out)**2.0) + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0)/(cos(i)**2.0 + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0)')

                
                #pars.add('y_tilt2', expr='-(180.0/pi)*asin((-cos(i)*sin(epsilon_in) + sin(i)*cos((z_tilt - z_tilt2)*(pi/180.0))*sqrt((cos(i)**2.0 - sin(epsilon_in)**2.0) + sin(i)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0))/(cos(i)**2.0 + sin(i)**2.0*cos((z_tilt-z_tilt2)*(pi/180.0))**2.0))')



                
                #pars.add('epsilon_out', True, epsilon_out_min, epsilon_out_max, expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt*(pi/180.0))*cos(z_tilt*(pi/180.0)))')
                #pars.add('epsilon_in', True, epsilon_in_min, epsilon_in_max, expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt2*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt2*(pi/180.0))*cos(z_tilt*(pi/180.0) - z_tilt2*(pi/180.0)))')

                #pars.add('y_tilt', expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt*(pi/180.0))*cos(z_tilt*(pi/180.0))) >= epsilon_out')
                #pars.add('y_tilt2', expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt2*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt2*(pi/180.0))*cos(z_tilt*(pi/180.0) - z_tilt2*(pi/180.0))) >= epsilon_in')

                
                pars.add('disk_flux', expr='0')







            def residual(pars, x, data, weight):

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
                #print '{:8f}'.format(np.sum(np.square(data_short - f(x_short)))/N), '\t', '{:8f}'.format(z_tilt), '\t', '{:8f}'.format(y_tilt),  '\t','{:8f}'.format(y_tilt2), '\t', '{:8f}'.format(z_tilt2), '\t', '{:10f}'.format(T_disk), '\t', '{:E}'.format(Lx)
                #print np.sum(np.square(data_short - f(x_short)))/N
        
                return (data_short - f(x_short))/weight








            

            mini = lmfit.Minimizer(residual, pars, fcn_args=(x, data, weight))
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
