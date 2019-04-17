# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import *
rcParams['mathtext.default'] = 'regular'
import math
from scipy import signal
import scipy

def main():	
	
	times = np.linspace(0,256e-9,256)
	b, a = signal.butter(9, 0.6)
	
	noise_1 = np.random.normal(0,11*1e-6,256)
	noise_2 = np.random.normal(0,11*1e-6,256)

	imp_1 = signal.unit_impulse(256, 100)
	imp_2 = signal.unit_impulse(256, 150)

	response_1 = signal.lfilter(b, a, imp_1)
	response_2 = signal.lfilter(b, a, imp_2)

	noisy_signal_1 = (response_1*4*1e-6)+noise_1
	noisy_signal_2 = (response_2*4*1e-6)+noise_2
	
	fig = plt.figure(figsize=(11,8.5)) #make a figure object
	ax = fig.add_subplot(2,1,1) #make a subplot for the limit
	ax.plot(times,noisy_signal_1)
	ax = fig.add_subplot(2,1,2) #make a subplot for the limit
	ax.plot(times,noisy_signal_2)
	fig.savefig("test.png",edgecolor='none',bbox_inches="tight") #save the figure
		
main()