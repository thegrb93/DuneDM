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
dir_state = [1, 1, 1, 1]

index = 0
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
	inputFilename = 'results/chimass_%.3f_%.3f_%.3f_%.3f.root' % (vpmass, chimass, kappa, alpha)
	subprocess.call(['cp', 'Events/chimass/unweighted_events.root', 'results/'+outputFilename])
	i = 0
	while i<4:
		iteration_state[i] += dir_state[i]
		if iteration_state[i]>=len(param_ranges[i]) or iteration_state[i]<0:
			dir_state[i] *= -1
			iteration_state[i] = 0
			i += 1
		else:
			break