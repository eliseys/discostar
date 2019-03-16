from __future__ import print_function

import lmfit
import numpy as np
import discostar_5_lib as d5l
import datetime
import itertools
import json
from collections import OrderedDict
import subprocess


data = json.load(open('data.json'), object_pairs_hook=OrderedDict)

p = json.load(open('parameters.json'), object_pairs_hook=OrderedDict)

iterables = tuple(p['additional'].values()) + tuple(p['main'].values())

output_file = './logs/' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S") + str() + ".json"

output = [p]

#for n in data:
for n in ['0.0']:
    for t in itertools.product(*iterables):

        parameters = dict(zip(tuple(p['additional'].keys()) + tuple(p['main'].keys()), t))
        parameters['n_data'] = n
        
        parameters = d5l.z_tilt_model(parameters)
        
        lmfit_parameters = lmfit.Parameters()
        
                
        for x in p['lmfit']['parameters']:
            lmfit_parameters.add_many(p['lmfit']['parameters'][x])
        
            
        mini = lmfit.Minimizer(d5l.residual, lmfit_parameters, fcn_args=(parameters, data))

        result = mini.minimize(method = p['lmfit']['method'], epsfcn = 5.0E-5)

        curve = {}
        
        for parameter in lmfit_parameters.keys():
            parameters[parameter] = result.params[parameter].value

        parameters = d5l.z_tilt_model(parameters)
        
            
        for UBV_filter in data[n].keys():
            if data[n][UBV_filter]['flux'] == []:
                continue
        
            parameters['UBV_filter'] = UBV_filter
            f = d5l.lc(parameters)
            
            curve[UBV_filter] = {"orbit": tuple(np.linspace(0, 1, parameters['lc_num'])), "flux": tuple(f(np.linspace(0, 1, parameters['lc_num'])))}

        #print(curve)
        #print(parameters)

        parameters['epsilon_out'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(result.params['y_tilt'].value/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(result.params['y_tilt'].value/180.0))*np.cos(np.pi*(parameters['z_tilt']/180.0))) 

        parameters['epsilon_in'] = (180.0/np.pi)*np.arcsin(np.cos(np.pi*(parameters['inclination']/180.0))*np.cos(np.pi*(result.params['y_tilt2'].value/180.0)) + np.sin(np.pi*(parameters['inclination']/180.0))*np.sin(np.pi*(result.params['y_tilt2'].value/180.0))*np.cos(np.pi*(parameters['z_tilt2']/180.0))) 
        
        
        output_i = {'parameters':parameters, 'data':data[n], 'curve':curve, 'rchisq':result.redchi, 'lmfit_parameters':p['lmfit']['parameters'], 'lmfit_method':p['lmfit']['method']}


        output.append(output_i)


        
with open(output_file, "w") as f:
    json.dump(output, f)


