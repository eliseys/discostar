import lmfit
import numpy as np
import discostar4_lib as d4l

h = 0.023
R = 0.3

#wobbling constant


half_hr_deg = (0.5*h/R)*(180.0/np.pi)

phase = {
    '00_05.dat':[324,	 0,     15.0,  21.0, 1],
    '05_10.dat':[342,	 18,    15.0,  21.0, 1],
    '10_15.dat':[0,	 36,    15.0,  21.0, 1],
    '15_20.dat':[18,     54,    15.0,  21.0, 1],
    '20_25.dat':[36,	 72,    15.0,  21.0, 1],
    '25_30.dat':[54,	 90,    15.0,  30.0,  1],
    '30_35.dat':[72,	 108,   15.0,  30.0,  1],
    '35_40.dat':[90,	 126,   15.0,  30.0,  1],
    '40_45.dat':[108,	 144,   15.0,  20.0, 0],
    '45_50.dat':[126,	 162,   10.0,  20.0, -1],
    '50_55.dat':[144,	 180,   15.0,  16.0, -1],
    '55_60.dat':[162,	 198,   14.0,  15.0, -1],
    '60_65.dat':[180,	 216,   13.0,  14.0, -1],
    '65_70.dat':[198,	 234,   12.0,  14.0, -1],
    '70_75.dat':[216,	 252,   10.0,  12.0, -1],
    '75_80.dat':[234,	 270,   11.0,  13.0, -1],
    '80_85.dat':[252,	 288,   15.0,  16.0, -1],
    '85_90.dat':[270,	 306,   16.0,  17.0, 0],
    '90_95.dat':[288,	 324,   17.0,  18.0, 1],
    '95_00.dat':[306,	 342,   17.0,  21.0, 1],
}

theta_list = [-3.0,]
kappa_list = [15.0,]
#inclination_list = [82.0,] 

inclination = 87.0

inc = inclination*(np.pi/180.0)


for theta in theta_list:
    for kappa in kappa_list:
        for i, x in enumerate(phase):

            z_tilt = phase[x][0]
            
            PSI_pr = phase[x][1]
            
            #if PSI_pr < (342.0 - 36.0):
            #    PSI_pr = PSI_pr + 36.0
            #elif PSI_pr >= (342.0 - 36.0):
            #    PSI_pr = PSI_pr + 36.0 - 360.0
            
            #y_tilt_min = phase[x][2]
            #y_tilt_max = phase[x][3]
            #eps_in_index = phase[x][4]

            data_file = x
            
            print  z_tilt, PSI_pr, data_file
            
            path_to_data_file = './OBSERVED_DATA/'+data_file

            pars = lmfit.Parameters()

            with open(path_to_data_file, 'r') as data:
                x, data = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)
                
                pars.add_many(
                    ('border', 0.13, False),
                    ('Lx', 1.5e+37, True, 0.5e+37, 5.0e+37, None, 1.0e+37),
                    ('PSI_pr', PSI_pr, False, 0, 360, None),
                    ('z_tilt', z_tilt, False, 0, 360, None),
                    ('y_tilt', 15.0, True, 0.0, 30.0, None),
                    #('y_tilt2', 10.0, True, 0.0, 30.0, None, 2.5),
                    ('delta', 0.5, True, 0.0, 1.0, None),
                    #('z_tilt2', 40.0, True, z_tilt2_min, z_tilt2_max, None),
                    ('z_tilt2', -40.0, True, -110.0, 0.0, None),
                    ('T_disk', 15000.0, True, 10000.0, 40000.0, None),
                    #('inclination', inclination, False, 0, 90),
                    #('eps_in_index', eps_in_index, False)
                )
                
                #pars.add('y_tilt', expr='20.0 * (0.6+0.4*cos(((z_tilt - 270.0)/360.0) * 2.0 * pi + 15.0 * pi/180.0 - 0.25 * 2.0 * pi))')
                pars.add('y_tilt2', expr='y_tilt*delta')

                #pars.add('y_tilt2', expr='y_tilt')

                #pars.add('aux1', expr='y_tilt')
                #pars.add('aux2', expr='y_tilt')

                #pars.add('z_tilt2', expr='(180.0/pi)*acos(-(tan(inclination*(pi/180.0)))**(-1.0)*(tan(y_tilt2*(pi/180.0)))**(-1.0)) - z_tilt'),

            

                #pars.add('cos_theta_minus_out', expr='(cos(inc)*sin(epsilon_out*(pi/180.0)) - sin(inc)*cos(z_tilt*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_out*(pi/180.0))**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0)')
                #pars.add('cos_theta_plus_out', expr='(cos(inc)*sin(epsilon_out*(pi/180.0)) + sin(inc)*cos(z_tilt*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_out*(pi/180.0))**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos(z_tilt*(pi/180.0))**2.0)')
                #pars.add('cos_theta_minus_in', expr='(cos(inc)*sin(epsilon_in*(pi/180.0)) - sin(inc)*cos((z_tilt - z_tilt2)*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_in*(pi/180.0))**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0)')
                #pars.add('cos_theta_plus_in', expr='(cos(inc)*sin(epsilon_in*(pi/180.0)) + sin(inc)*cos((z_tilt - z_tilt2)*(pi/180.0)) * sqrt( cos(inc)**2.0 - sin(epsilon_in*(pi/180.0))**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0 ) )/(cos(inc)**2.0 + sin(inc)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0)')

                #pars.add('y_tilt', expr='min((180.0/pi)*acos(cos_theta_minus_out), (180.0/pi)*acos(cos_theta_plus_out))')

                #pars.add('y_tilt2', expr='min((180.0/pi)*acos(cos_theta_minus_in), (180.0/pi)*acos(cos_theta_plus_in))')

                
                #pars.add('y_tilt', expr='-(180.0/pi)*asin((-cos(i)*sin(epsilon_out) + sin(i)*cos(z_tilt*(pi/180.0))*sqrt((cos(i)**2.0 - sin(epsilon_out)**2.0) + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0))/(cos(i)**2.0 + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0))')
                #pars.add('y_tilt', expr='sqrt((cos(i)**2.0 - sin(epsilon_out)**2.0) + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0)/(cos(i)**2.0 + sin(i)**2.0*cos(z_tilt*(pi/180.0))**2.0)')

                
                #pars.add('y_tilt2', expr='-(180.0/pi)*asin((-cos(i)*sin(epsilon_in) + sin(i)*cos((z_tilt - z_tilt2)*(pi/180.0))*sqrt((cos(i)**2.0 - sin(epsilon_in)**2.0) + sin(i)**2.0*cos((z_tilt - z_tilt2)*(pi/180.0))**2.0))/(cos(i)**2.0 + sin(i)**2.0*cos((z_tilt-z_tilt2)*(pi/180.0))**2.0))')



                
                #pars.add('epsilon_out', False, expr='(180.0/pi)*asin(cos(inclination*(pi/180.0))*cos(y_tilt*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt*(pi/180.0))*cos(z_tilt*(pi/180.0)))')
                #pars.add('epsilon_in', False, expr='(180.0/pi)*asin(cos(inclination*(pi/180.0))*cos(y_tilt2*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt2*(pi/180.0))*cos(z_tilt*(pi/180.0) + z_tilt2*(pi/180.0)))')

                #pars.add('y_tilt', expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt*(pi/180.0))*cos(z_tilt*(pi/180.0))) >= epsilon_out')
                #pars.add('y_tilt2', expr='asin(cos(inclination*(pi/180.0))*cos(y_tilt2*(pi/180.0)) + sin(inclination*(pi/180.0))*sin(y_tilt2*(pi/180.0))*cos(z_tilt*(pi/180.0) - z_tilt2*(pi/180.0))) >= epsilon_in')

                
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

                #border = pars['border'].value

                x_short = []
                data_short = []

                border = 0.13

                #border2 = 0.31

                #border3 = 0.233
                #border4 = 0.33
                
                for i, X in enumerate(x):
                    if (X > border) or (X < 1.0 - border):
                    #if ((X > border) and (X < border1)) or ((X > border2) and (X < 1.0 - border)):
                        x_short.append(x[i])
                        data_short.append(data[i])

                data_short = np.array(data_short)

                N = len(data_short)
        
                f = d4l.lc_disk(z_tilt, PSI_pr, Lx, y_tilt, y_tilt2, z_tilt2, disk_flux, T_disk, theta, kappa)

                #if np.sum(np.square(data_short - f(x_short)))/N < 0.07:
                #
                #print '{:8f}'.format(np.sum(np.square(data_short - f(x_short)))/N), '\t', '{:8f}'.format(z_tilt), '\t', '{:8f}'.format(y_tilt),  '\t','{:8f}'.format(y_tilt2), '\t', '{:8f}'.format(z_tilt2), '\t', '{:10f}'.format(T_disk), '\t', '{:E}'.format(Lx)
                #print np.sum(np.square(data_short - f(x_short)))/N
        
                return data_short - f(x_short)





            
            mini = lmfit.Minimizer(residual, pars, fcn_args=(x, data))
            result = mini.minimize(epsfcn=1.0e-4)
            #result = mini.minimize(method='differential_evolution')
            #result = mini.minimize()


            #mini = lmfit.Minimizer(test_fit, pars, fcn_args=(x, data))
            #result = mini.minimize()


            #result = mini.minimize(method='differential_evolution')
            
            #print result.params
            
            #print result.params.pretty_print

    
            print(lmfit.fit_report(result))
            print "kappa", "\t", kappa
            print "ns_theta", "\t", theta
            #print(fit_report(mini))
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
