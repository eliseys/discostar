from __future__ import print_function

import lmfit
import numpy as np
import discostar_6_lib as d6l
import datetime
import itertools
import json
from collections import OrderedDict
import subprocess
from scipy.integrate import quad
from scipy import interpolate
import sys


def dot_psi(x, z_tilt_dot_profile, break_amplitude, break_duration, slowest_phase):
    if z_tilt_dot_profile == 1:
        return 1.0 - break_amplitude/(2.0*np.pi) * np.exp((-(((x - slowest_phase + 0.5)%1.0-0.5))**2.0)/(2.0*break_duration**2.0))
    elif z_tilt_dot_profile == 2:
        return 1.0 - break_amplitude * np.sin((x - slowest_phase + 0.25)*2*np.pi)
    elif z_tilt_dot_profile == 3:
        b_slope = 2.5
        return 1.0 - break_amplitude/(1.0 + np.abs(((x - slowest_phase + 0.5)%1.0-0.5)/break_duration)**(2.0*b_slope))


def psi(x, z_tilt_dot_profile, break_amplitude, break_duration, slowest_phase, max_opening_phase):
    normalizaton_coeff = quad(dot_psi, 0, 1, args=(z_tilt_dot_profile, break_amplitude, break_duration, slowest_phase))[0]
    max_opening_phase_psi = quad(dot_psi, 0, max_opening_phase, args=(z_tilt_dot_profile, break_amplitude, break_duration, slowest_phase))[0]

    return 360.0 * (quad(dot_psi, 0, x, args=(z_tilt_dot_profile, break_amplitude, break_duration, slowest_phase))[0]/normalizaton_coeff - max_opening_phase_psi)



data = json.load(open('DATA/data.json'), object_pairs_hook=OrderedDict)
N = 20.0 # precession cycle lenght in orbital periods

p = json.load(open('parameters.json'), object_pairs_hook=OrderedDict)

iterables = tuple(p['additional'].values()) + tuple(p['main'].values())

output_file = './logs/' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S") + str() + ".json"

output = [p]


for n in data:
#for n in ['0']:

    #if len(data[n]['B']['MJD']) + len(data[n]['V']['MJD']) < 15.0: #minimal number of points at lightcurve
    #   continue
    if len(data[n]['B']['JD']) + len(data[n]['V']['JD']) < 15.0: #minimal number of points at lightcurve
        continue

    for t in itertools.product(*iterables):


        parameters = dict(zip(tuple(p['additional'].keys()) + tuple(p['main'].keys()), t))
        parameters['n_data'] = n
        parameters['number_of_orbits'] = N
        
        parameters = d6l.z_tilt_model(parameters)

        lmfit_parameters = lmfit.Parameters()

        #parameters['inclination'], parameters['y_tilt'] = (360.0/(2.0*np.pi)) * d6l.i_and_theta_calculator(parameters['max_opening_phase'], np.arctan(0.5*parameters['h']/parameters['R']))
        parameters['inclination'], y_tilt_initial = (360.0/(2.0*np.pi)) * d6l.i_and_theta_calculator(parameters['max_opening_phase'], np.arctan(0.5*parameters['h']/parameters['R']))

        #print(parameters['inclination'], y_tilt_initial, parameters['h'])
        #sys.exit()

        
        if parameters['inclination'] < 85.0:
            continue
        
        a_list = [6.525e11, 6.516168e11, 6.509544e11, 6.504936e11, 6.502152e11, 6.501e11]
        i_list = [85.0, 86.0, 87.0, 88.0, 89, 90.0]
        
        f = interpolate.interp1d(i_list, a_list, kind='cubic')
        
        parameters['a'] = float(f(parameters['inclination']))

        
        for x in p['lmfit']['parameters']:
            lmfit_parameters.add_many(p['lmfit']['parameters'][x])

        #z_tilt_average = psi((1.0/len(data))*float(n), parameters["z_tilt_dot_profile"], parameters["break_amplitude"], parameters["break_duration"], parameters["slowest_phase"], parameters["max_opening_phase"])


        a = 360.0/N
        b = - a * parameters['max_opening_phase'] * N
        z_tilt_linear = a * int(float(n)) + b # degrees
        
        y_tilt = (u'y_tilt', y_tilt_initial, True, 5.0, 25.0, None, None)
        y_tilt2 = (u'y_tilt2', y_tilt_initial, True, 0.0, 25.0, None, None)

        z_tilt = (u'z_tilt', z_tilt_linear, True, z_tilt_linear - 20.0, z_tilt_linear + 20.0, None, None)
        
        lmfit_parameters.add_many(y_tilt)
        lmfit_parameters.add_many(y_tilt2)
        lmfit_parameters.add_many(z_tilt)


        #lmfit_parameters.add_many(y_tilt)
        

        #print(lmfit_parameters)

        #sys.exit()
        
        mini = lmfit.Minimizer(d6l.residual, lmfit_parameters, fcn_args=(parameters, data), nan_policy = 'omit')
        
        result = mini.minimize(method = p['lmfit']['method'], epsfcn = 5.0E-5, max_nfev = 0)
        #result = mini.minimize(method = p['lmfit']['method'])

        curve = {}
       
        for parameter in lmfit_parameters.keys():
            parameters[parameter] = result.params[parameter].value

        parameters = d6l.z_tilt_model(parameters)

            
        for UBV_filter in data[n].keys():
            if data[n][UBV_filter]['flux'] == []:
                continue
        
            parameters['filter'] = UBV_filter
            parameters['lc_start'] = 0.13
            parameters['lc_end'] = 0.87

            x, y = d6l.lc(parameters)
            
            curve[UBV_filter] = {"orbit": list(x), "flux": list(y)}


        # parameters['epsilon_out'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(result.params['y_tilt'].value/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(result.params['y_tilt'].value/180.0))*np.cos(np.pi*(parameters['z_tilt']/180.0))) 

        # parameters['epsilon_in'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(result.params['y_tilt2'].value/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(result.params['y_tilt2'].value/180.0))*np.cos(np.pi*(parameters['z_tilt2']/180.0))) 

        parameters['epsilon_out'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(parameters['y_tilt']/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(parameters['y_tilt']/180.0))*np.cos(np.pi*(parameters['z_tilt']/180.0))) 

        parameters['epsilon_in'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(parameters['y_tilt2']/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(parameters['y_tilt2']/180.0))*np.cos(np.pi*(parameters['z_tilt2']/180.0))) 

        #parameters['epsilon_out'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(parameters['y_tilt']/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(parameters['y_tilt']/180.0))*np.cos(np.pi*(parameters['z_tilt']/180.0))) 

        #parameters['epsilon_in'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(result.params['y_tilt2'].value/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(result.params['y_tilt2'].value/180.0))*np.cos(np.pi*(parameters['z_tilt2']/180.0))) 
        
        output_i = {'parameters':parameters, 'data':data[n], 'curve':curve, 'rchisq':result.redchi, 'lmfit_parameters':p['lmfit']['parameters'], 'lmfit_method':p['lmfit']['method']}
        

        output.append(output_i)
        
data = json.load(open('DATA/data.json'), object_pairs_hook=OrderedDict)
json.dump(output, f)


