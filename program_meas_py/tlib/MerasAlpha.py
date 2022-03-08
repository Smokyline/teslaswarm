import numpy as np
from program_meas_py.tlib.cmp import *

gamma = 1.0
p_cmp = n_cmp1

def n_cmp_alpha(zn, level, alpha):

	tmp = p_cmp(zn, level) - alpha;
	if (tmp >= 0.0):
		tmp /= (1-alpha)
	else:
		tmp /= (1+alpha)

	return tmp

def MerasAlpha(l_len, weight, vscale, value, talpha=None):
	
	if not talpha is None:
		if (talpha > 0):
			n_fun = int(talpha)
			alpha = talpha - n_fun
		else:
			n_fun = int(abs(talpha))
			alpha = talpha + n_fun

		if n_fun == 1:
			p_cmp = n_cmp1
		elif n_fun == 2:
			p_cmp = n_cmp2
		elif n_fun == 3:
			p_cmp = n_cmp3
		elif n_fun == 4:
			p_cmp = n_cmp4
		elif n_fun == 5:
			p_cmp = n_cmp5
		elif n_fun == 6:
			p_cmp = n_cmp6

	mera = np.zeros_like(value)

	if alpha == 0:
		for i in range(len(weight)):
			if (weight[i] == 0): continue
			for j in range(len(value)):
#				n_cmp_alpha(vscale[i], value, alpha)
				mera[j] += n_cmp_alpha(vscale[i], value[j], alpha)*weight[i]
	else:
		for i in range(len(weight)):
			if (weight[i] == 0): continue
			for j in range(len(value)):
				mera[j] += n_cmp_alpha(vscale[i], value[j], alpha)*weight[i]

	mera /= l_len

	return mera