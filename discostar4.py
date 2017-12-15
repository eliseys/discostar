import lmfit
import numpy as np
import discostar_4_lib as d4l

data_file = '75_80.dat'
path_to_data_file = './OBSERVED_DATA/'+data_file

pars_disk = lmfit.Parameters()
pars = lmfit.Parameters()

with open(path_to_data_file, 'r') as data:
    x, data = np.loadtxt(data, dtype=('float'), usecols=(4,1), unpack=True)

pars_disk.add_many(
    ('border', 0.15, False),
    ('Lx', 6.46e+37, False),
    ('PSI_pr', 90, False),
    ('z_tilt', 180.0, False),
    ('y_tilt', 15.0, False),
    ('y_tilt2', 13.4, False),
    ('z_tilt2', 0.0, False),
    ('T_disk', 17000.0, True, 10000.0, 20000.0))


def residual_disk(pars_disk, x, data):

    z_tilt = pars_disk['z_tilt'].value
    PSI_pr = pars_disk['PSI_pr'].value

    y_tilt = pars_disk['y_tilt'].value
    T_disk = pars_disk['T_disk'].value
    border = pars_disk['border'].value

    x_short = []
    data_short = []
    
    for i, X in enumerate(x):
        if (X < border) or (X > 1.0 - border):
            x_short.append(x[i])
            data_short.append(data[i])
        
    f = d4l.lc_disk(y_tilt, T_disk, z_tilt, PSI_pr)

    return data_short - f(x_short)


min_disk = lmfit.Minimizer(residual_disk, pars_disk, fcn_args=(x, data))
result_disk = min_disk.minimize()

print result_disk.params


#pars.add_many(
#    ('border', 0.15, False),
#    ('Lx', 6.46e+37, True, 1.0e+37, 9.0e+37),
#    ('PSI_pr', 90, False),
#    ('z_tilt', 180.0, False),
#    ('y_tilt', result_disk.params['y_tilt'], True, 0.0, 30.0),
#    ('y_tilt2', 13.4, True, 0.0, 30.0),
#    ('z_tilt2', 0.0, True, 0.0, 90.0),
#    ('T_disk', result_disk.params['T_disk'], True, 10000.0, 20000.0))




pars.add_many(
    ('border', 0.15, False),
    ('Lx', 6.46e+37, False, 1.0e+37, 9.0e+37),
    ('PSI_pr', 90, False),
    ('z_tilt', 180.0, False),
    ('y_tilt', result_disk.params['y_tilt'].value, True, 0.0, 30.0),
    ('y_tilt2', 13.4, False, 0.0, 30.0),
    ('z_tilt2', 0.0, False, 0.0, 90.0),
    ('T_disk', result_disk.params['T_disk'].value, True, 10000.0, 20000.0))


def residual(pars, x, data):

    z_tilt = pars['z_tilt'].value
    PSI_pr = pars['PSI_pr'].value
    
    Lx = pars['Lx'].value
    y_tilt2 = pars['y_tilt2'].value
    z_tilt2 = pars['z_tilt2'].value

    y_tilt = pars['y_tilt'].value
    T_disk = pars['T_disk'].value
    
    f = d4l.lc(Lx, y_tilt2, z_tilt2, y_tilt, T_disk, z_tilt, PSI_pr)

    return data - f(x)


min_ = lmfit.Minimizer(residual, pars, fcn_args=(x, data))
result = min_.minimize()



print '\n'


#print result.params.pretty_print

print(lmfit.fit_report(result.params))

ci = lmfit.conf_interval(min_, result)
lmfit.printfuncs.report_ci(ci)
