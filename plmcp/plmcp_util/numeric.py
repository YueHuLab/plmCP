'''numerical array calculations powered by numba'''

from typing import Union, List, Tuple, Dict

import numpy as np
import numba

#@numba.njit('f4[:,:],f4[:,:](f4[:,:], f4)', nogil=True, fastmath=True, cache=True)
@numba.njit( nogil=True, fastmath=True, cache=True)
def fill_matrix_local(a: np.ndarray, gap_extension : np.float32 ):
	nrows: int = a.shape[0] + 1
	ncols: int = a.shape[1] + 1
	H: np.ndarray = np.zeros((nrows, ncols), dtype=np.float32)
	Antecedent: np.ndarray = np.zeros((nrows, ncols), dtype=np.float32)
	h_tmp: np.ndarray = np.zeros(4, dtype=np.float32)
	for i in range(1, nrows):
		for j in range(1, ncols):
			h_tmp[0] = H[i-1, j-1] + a[i-1, j-1]
			h_tmp[1] = H[i-1, j] - gap_extension
			h_tmp[2] = H[i, j-1] - gap_extension
			H[i, j] = np.max(h_tmp)
			if (H[i,j]==h_tmp[0]):
                            Antecedent[i,j]=1.0
			if (H[i,j]==h_tmp[1]):
                            Antecedent[i,j]=2.0
			if (H[i,j]==h_tmp[2]):
                            Antecedent[i,j]=3.0
			#print("HH 0 \n")
			#print(h_tmp[0])
			#print("HH 1 \n")
			#print(h_tmp[1])
			#print("HH 2 \n")
			#print(h_tmp[2])
			#print("H \n")
			#print(H[i, j])
			#print("\n")
	#print("H \n")
	#print(H)
	#print("H \n")
	#return H, Antecedent
	return H, Antecedent


@numba.njit('f4[:,:](f4[:,:], f4)', nogil=True, fastmath=True, cache=True)
def fill_matrix_global(a: np.ndarray, gap_extension : np.float32):
	'''
	fill score matrix in Needleman-Wunch procedure - global alignment
	Params:
		a: (np.array)
		gap_penalty (float)
	Return:
		b: (np.array)
	'''
	nrows: int = a.shape[0] + 1
	ncols: int = a.shape[1] + 1
	H: np.ndarray = np.zeros((nrows, ncols), dtype=np.float32)
	h_tmp: np.ndarray = np.zeros(3, dtype=np.float32)
	for i in range(0, nrows):
		for j in range(0, ncols):
			if ((i==0) and (j==0)):
				H[i, j] = 0
			elif ((i==0) or (j==0)):
				H[i, j] = - (i+j-1) * gap_extension
			else:
				h_tmp[0] = H[i-1, j-1] + a[i-1, j-1]
				h_tmp[1] = H[i-1, j] - gap_extension
				h_tmp[2] = H[i, j-1] - gap_extension
				H[i, j] = np.max(h_tmp)
	#print("H \n")
	#print(H)
	#print("H \n")

	return H


def fill_score_matrix(sub_matrix: np.ndarray,
					  gap_extension: np.float32, #huyue
					  mode: str = 'local') -> np.ndarray:
	'''
	use substitution matrix to create score matrix
	set mode = local for Smith-Waterman like procedure (many local alignments)
	and mode = global for Needleamn-Wunsch like procedure (one global alignment)
	Params:
		sub_matrix: (np.array) substitution matrix in form of 2d
			array with shape: [num_res1, num_res2]
		gap_penalty: (float)
		mode: (str) set global or local alignment procedure
	Return:
		score_matrix: (np.array)
	'''
	assert gap_extension >= 0, 'gap extension must be positive'
	assert isinstance(mode, str)
	assert mode in {"global", "local"}
	assert isinstance(gap_extension, (int, float))
	assert isinstance(sub_matrix, np.ndarray), \
		'substitution matrix must be numpy array'
	# func fill_matrix require np.float32 array as input
	if not np.issubsctype(sub_matrix, np.float32):
		sub_matrix = sub_matrix.astype(np.float32)
	if mode == 'local':
		score_matrix,antecedent = fill_matrix_local(sub_matrix, gap_extension = gap_extension)
	elif mode == 'global':
		score_matrix = fill_matrix_global(sub_matrix, gap_extension = gap_extension)
	return score_matrix,antecedent


@numba.njit('types.Tuple((f4, i4))(f4, f4, f4)', cache=True)
def max_from_3(x: float, y: float, z: float) -> Tuple[float, int]:
	'''
	return value and index of biggest values
	'''
	# 2 idx should be diagonal
	if z >= y and z >= x:
		return z, 2
	if x > y and x > z:
		return x, 0
	else:
		return y, 1


@numba.jit(fastmath=True, cache=True)
def traceback_from_point_opt2(scoremx: np.ndarray, antecedent: np.ndarray, point: Tuple[int, int],gap_extension: np.float32, #huyue#huyue 1e-3
							mode: str = 'local', stop_value: np.float32 = 1e-5) -> np.ndarray:
	'''
	find optimal route over single path
	Args:
		scoremx (np.ndarray 2D):
		point (tuple): y, x coordinates
		stop_value (float): end of route criteria
	Returns:
		ndarray coordinates of path
	'''
	assert isinstance(mode, str)
	assert mode in {"global", "local"}
	f_right: np.float32 = 0.0
	f_left: np.float32 = 0.0
	#f_diag: int = 0
	#fi_max: int = 0
	f_diag: np.float32 = 0.0 #huyue
	fi_max: np.float32 = 0.0
	fi_max_old: np.float32 = 0.0
	fi_score: np.float32 = 0.0
	testValue: np.float32 = 0.0
	# assume that the first move through alignment is diagonal
	fi_argmax: int = 2
	fi_argmax_old: int = 2
	y_size: int = scoremx.shape[0]
	x_size: int = scoremx.shape[1]
	yi: int = point[0]
	xi: int = point[1]
	assert y_size > yi
	assert x_size > xi
	# set starting position
	position: int = 1
	# maximum size of path
	size: int = y_size + x_size
	path_arr: np.ndarray = np.zeros((size, 2), dtype=np.int32)
	# do not insert starting point
	path_arr[0, 0] = yi
	path_arr[0, 1] = xi
	# iterate until border is hit
	# score matrix have one extra row and column
	while ((yi > 1) or (xi > 1)):
		# find previous fi_argmax was diagnal
		#if (abs(antecedent[yi,xi]-0.0)<stop_value):
		if (antecedent[yi,xi]==0.0):
                    break
		if (antecedent[yi,xi]==1.0):
                    fi_max=scoremx[yi-1, xi-1]
                    fi_argmax=2
		if (antecedent[yi,xi]==2.0):
                    fi_max=scoremx[yi-1, xi]
                    fi_argmax=0
		if (antecedent[yi,xi]==3.0):
                    fi_max=scoremx[yi, xi-1]
                    fi_argmax=1
		fi_score=scoremx[yi,xi]

                # add point to path
			# diagonal move
		if fi_argmax == 2:
			yi_new = yi - 1
			xi_new = xi - 1
		# move left
		elif fi_argmax == 1:
			yi_new = yi
			xi_new = xi - 1
		# move right
		else:
			yi_new = yi - 1
			xi_new = xi
		# store index
		path_arr[position, 0] = yi_new
		path_arr[position, 1] = xi_new
		#huyue
                # set new indices
		yi = yi_new
		xi = xi_new
		position += 1
		# if maximal value if <= 0 stop loop
		if ((mode == 'local') and (fi_max <= stop_value)):
			break
		
                #if ((mode == 'local') and (fi_argmax==0 and (abs(fi_score-fi_max)-gap_extension))):
		#	break
		if ((mode == 'local') and (fi_score <= stop_value)):
			break
	        
	# push one index up to remove zero padding effect
	# not done
	#print(f"path_arr\t {path_arr}\n")
	#path_arr = path_arr[:position, :]
	path_arr = path_arr[:position-1, :]
	return path_arr
