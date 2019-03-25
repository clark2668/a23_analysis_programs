# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import *
rcParams['mathtext.default'] = 'regular'
import math
from scipy import signal
import scipy

td_args = {
        'down1': (-0.8, 15e-9, 2.3e-9, 0),
        'down2': (-0.2, 15e-9, 4e-9, 0),
        'up': (1, 18e-9, 7e-9, 1e9)
}

td_args['up'] = (-np.sqrt(2 * np.pi) * (td_args['down1'][0] * td_args['down1'][2] +
                       td_args['down2'][0] * td_args['down2'][2]) /
					(2e18 * td_args['up'][2] ** 3),) + td_args['up'][1:]
					
def td_fdown1(x):
        return (td_args['down1'][3] + td_args['down1'][0] *
                np.exp(-(x - td_args['down1'][1]) ** 2 /
			(2 * td_args['down1'][2] ** 2)))
			
def td_fdown2(x):
        return (td_args['down2'][3] + td_args['down2'][0] *
                np.exp(-(x - td_args['down2'][1]) ** 2 /
			(2 * td_args['down2'][2] ** 2)))
			
def td_fup(x):
        return (td_args['up'][0] *
                (td_args['up'][3] * (x - td_args['up'][1])) ** 2 *
			np.exp(-(x - td_args['up'][1]) / td_args['up'][2]))


def tunnel_diode(self, dt, trace):
        t_max = 1e-7
        n_pts = int(t_max * dt)
        times = np.linspace(0, t_max, n_pts + 1)
        diode_resp = td_fdown1(times) + td_fdown2(times)
        t_slice = times > td_args['up'][1]
        diode_resp[t_slice] += td_fup(times[t_slice])
        conv = scipy.signal.convolve(trace ** 2 , diode_resp, mode='full')
        # conv multiplied by dt so that the amplitude stays constant for
        # varying dts (determined emperically, see ARVZAskaryanSignal comments)
        # Setting output
        trace_after_tunnel_diode = conv / dt
        trace_after_tunnel_diode = trace_after_tunnel_diode[:trace.shape[0]]
        return trace_after_tunnel_diode


def main():	
	
	times = np.linspace(0,256e-9,256)
	
	noise = np.random.normal(0,11*1e-6,256)
	imp = signal.unit_impulse(256, 'mid')
	b, a = signal.butter(9, 0.6)
	response = signal.lfilter(b, a, imp)
	noisy_signal = (response*300*1e-6)+noise
	td = tunnel_diode(noisy_signal,1e-9,noisy_signal)
	
	fig = plt.figure(figsize=(11,8.5)) #make a figure object
	ax = fig.add_subplot(2,1,1) #make a subplot for the limit
	ax.plot(times,noisy_signal)
	ax2 = fig.add_subplot(2,1,2) #make a subplot for the limit
	ax2.plot(times,td)
	fig.savefig("test.png",edgecolor='none',bbox_inches="tight") #save the figure
	
	times_new = np.linspace(0,256e-9,256)
	amp_levels = np.logspace(-3,0)
	#amp_levels = np.array([1e-7,1e-6,1e-7,1e-5,1e-3,1e-2,1e-1,1e0,1e1])
	max_levels = np.empty(len(amp_levels))
	volt_squared = np.empty(len(amp_levels))
	multiple = np.empty(len(amp_levels))
	for i, level in enumerate(amp_levels):
		max_levels[i] = np.max(np.abs(response*level))
		volt_squared[i]=level**2
		multiple[i]=max_levels[i]/volt_squared[i]
		
		
	print amp_levels[2]/amp_levels[0]	
	fig2 = plt.figure(figsize=(11,8.5)) #make a figure object
	ax3 = fig2.add_subplot(1,1,1) #make a subplot for the limit
	ax3.plot(amp_levels,max_levels,'o-')
	ax3.plot(amp_levels,volt_squared,'s-')
	#ax3.set_yscale('log')
	ax3.set_xscale('log')
	#ax3.axvline(x=11e-6,linestyle='--',color='red')
	fig2.savefig("test2.png",edgecolor='none',bbox_inches="tight") #save the figure
	
		
main()