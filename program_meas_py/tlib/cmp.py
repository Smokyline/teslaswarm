import numpy as np

def DOUBLE_EQ(A, B):
	return A==B

def n_cmp0(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0) and DOUBLE_EQ(B-A, 0)):
		return 0
	else:
		return (B - A)
	
def n_cmp1(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return (B - A) / np.maximum(A, B)
	
def n_cmp2(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return (B - A) / (np.abs(A) + np.abs(B - A))
	
def n_cmp3(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return (B - A) / (1.0 + np.abs(B - A))
	
def n_cmp4(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return 2.0 * np.arctan(B - A) / np.pi
	
def n_cmp5(A, B):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return (B - A) / (A + B)
	
def n_cmp6(A, B, gamma=1.0):
	if (DOUBLE_EQ(A, 0) and DOUBLE_EQ(B, 0)):
		return 0
	else:
		return (B - A) / np.power((np.power(A, 1.0/gamma) + np.power(B, 1.0/gamma)), gamma)
