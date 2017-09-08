import logging
import logging.config
import os
import re
import shutil
import subprocess
import sys
import time

param_ranges = [
[x * 0.1 for x in range(30, 80)], #vpmass
[x * 0.1 for x in range(1, 20)], #chimass
[x * 0.01 for x in range(1, 2)], #kappa
[x * 1 for x in range(1, 2)] #alpha
]
iteration_state = [0, 0, 0, 0]

param_card_f = open("Cards/param_card.dat", "r")
param_card = param_card_f.read()
param_card_f.close()

def setParams(vpmass, chimass, kappa, alpha):
	file_string = re.sub("\n.+?VPMASS\n", "\n        32     %.8e   # VPMASS\n"%vpmass, param_card)
	file_string = re.sub("\n.+?CHIMASS\n", "\n        33     %.8e   # CHIMASS\n"%chimass, file_string)
	
	file_handle = open("Cards/param_card.dat", 'w')
	file_handle.write(file_string)
	file_handle.close()

i = 0
while i<4:
	vpmass = param_ranges[0][iteration_state[0]]
	chimass = param_ranges[1][iteration_state[1]]
	kappa = param_ranges[2][iteration_state[2]]
	alpha = param_ranges[3][iteration_state[3]]
	print(vpmass, chimass, kappa, alpha)
	setParams(vpmass, chimass, kappa, alpha)
	subprocess.call('rm Events/chimass/*', shell = True)
	subprocess.call([sys.executable, '-O', 'bin/generate_events', '2', '2', 'chimass'])
	outputFilename = 'chimass_%.3f_%.3f_%.3f_%.3f.root' % (vpmass, chimass, kappa, alpha)
	subprocess.call(['cp', 'Events/chimass/unweighted_events.root', 'results/'+outputFilename])
	i = 0
	while i<4:
		iteration_state[i] += 1
		if iteration_state[i]>=len(param_ranges[i]):
			iteration_state[i] = 0
			i += 1
		else:
			break