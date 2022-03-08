import numpy as np

def CmpForm(level, weight):

	weight[:] = 0

	for i in range(len(level)):
		j = int(level[i]) - 1
		weight[j] = weight[j] + 1

	return weight