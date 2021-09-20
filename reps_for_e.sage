paramset_NTRU1 = {'n': 509, 'q': 2048, 'w': 254, 'dg':339}
paramset_NTRU2 = {'n': 677, 'q': 2048, 'w': 254, 'dg':451}
paramset_NTRU3 = {'n': 821, 'q': 4096, 'w': 510, 'dg':547}
paramset_NTRU_HRSS = {'n': 701, 'q': 8192, 'w': 468, 'dg':547}


paramset_NTRUPrime1 = {'n': 653, 'q': 4621, 'w': 288, 'dg':435}
paramset_NTRUPrime2 = {'n': 761, 'q': 4591, 'w': 286, 'dg':507}
paramset_NTRUPrime3 = {'n': 857, 'q': 5167, 'w': 322, 'dg':571}

paramset_BLISS = {'n': 512, 'q': 12289, 'w': 308, 'dg': 308}
paramset_GLP = {'n': 512, 'q': 8383489, 'w': 342, 'dg': 342}


NTRU_ieee401 = {'n': 401, 'q': 2048, 'w': 226, 'dg': 266}
NTRU_ieee449 = {'n': 449, 'q': 2048, 'w': 268, 'dg': 298}
NTRU_ieee677 = {'n': 677, 'q': 2048, 'w': 314, 'dg': 450}
NTRU_ieee1087 = {'n': 1087, 'q': 2048, 'w': 240, 'dg': 724}
NTRU_ieee541 = {'n': 541, 'q': 2048, 'w': 98, 'dg': 360}
NTRU_ieee613 = {'n': 613, 'q': 2048, 'w': 110, 'dg': 408}
NTRU_ieee887 = {'n': 887, 'q': 2048, 'w': 162, 'dg': 490}
NTRU_ieee1171 = {'n': 1171, 'q': 2048, 'w': 212, 'dg': 780}
NTRU_ieee659 = {'n': 659, 'q': 2048, 'w': 76, 'dg': 438}
NTRU_ieee761 = {'n': 761, 'q': 2048, 'w': 84, 'dg': 506}
NTRU_ieee1087_ = {'n': 1087, 'q': 2048, 'w': 126, 'dg': 724}
NTRU_ieee1499 = {'n': 1499, 'q': 2048, 'w': 158, 'dg': 998}


all_param_sets = [paramset_NTRU1, paramset_NTRU2, paramset_NTRU3,paramset_NTRU_HRSS,
				paramset_NTRUPrime1, paramset_NTRUPrime2, paramset_NTRUPrime3,
				 paramset_BLISS,paramset_GLP]
all_NTRU_ieee = [NTRU_ieee401, NTRU_ieee449, NTRU_ieee677, NTRU_ieee1087,
				NTRU_ieee541, NTRU_ieee613, NTRU_ieee887, NTRU_ieee1171,
				NTRU_ieee659, NTRU_ieee761, NTRU_ieee1087_, NTRU_ieee1499]


def multinom(n, c):

	"""
	Product of binimial coefficients:
	(n choose c[0]) * (n-c[0] choose c[1])*...*
	"""
	assert sum(c) == n, 'bad input to multinom!'
	res = 1
	n_ = n
	for i in range(len(c)):
		#print('from multinom:', n_, c[i], binomial(n_, c[i]))
		#if (i==len(c)-1):
		#	assert(binomial(n_, c[i])==1)
		res*=binomial(n_, c[i])
		n_ = n_ - c[i]
	return res


def search_space(param):
	n = param['n']
	w = param['w']
	return log(multinom(n,[round(w/2), round(w/2), n - 2*round(w/2)]),2).n()


def MitMError(params):
	n = params['n']
	q = params['q']
	w = params['w']

	min_runtime = float("infinity")
	Reps = (multinom(round(w/2), [round(w/4),round(w/2) - round(w/4) ]))**2
	k = round(log(Reps, q))
	L2_s = multinom(ceil(n/2), [ceil(w/8),ceil(w/8), ceil(n/2) - 2*ceil(w/8)])
	w_e = round(k*w/n) # expected weight of e on k coordinates
	assert(k/2>2*round(w_e/4))
	L2_e = multinom(k/2, [round(w_e/4),round(w_e/4), k/2 - 2*round(w_e/4)  ])

	L2 = L2_s*L2_e

	L1 = max(1,L2*L2 / Reps)
	L0 = max(1,L1**2 / 2**(floor(n-k)))

	runtime_list = [log(L2,2).n(), log(L1,2).n(), log(L0,2).n()]
	runtime = max(runtime_list)
	lvl = runtime_list.index(runtime)


	return runtime, lvl, runtime_list, k

"""
print("--------- MitM for error Rep-0---------")
for param in all_param_sets:
	res1 = MitMError(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', res1)
"""
# ------REP0------
def Rep0Error(params):
	n = params['n']
	q = params['q']
	w = params['w']
	we = params['dg']

	min_runtime = float("infinity")

	Reps_s = (multinom(round(w/2), [round(w/4),round(w/2) - round(w/4) ]))**2
	L2_original = multinom(ceil(n/2), [ceil(w/8),ceil(w/8), ceil(n/2) - 2*ceil(w/8)])
	k_s = round(log(Reps_s, q))

	#if k_s %2 == 1:
	#	k_s -=1

	#brute-force search for k (# of coordinates we merge on) that satisfies
	# k = k_s + log(Reps_e,2), where Reps_e is the number of representations
	# for the error on k coordinates
	err_in_solution = 100
	for k in range(k_s, round(n/2), 2):
		# expected error weight on k
		w_k = round(k*we/n)
		if w_k %2 == 1:
			w_k+=1

		R_e = max(1, (multinom(round(w_k/2), [round(w_k/4),round(w_k/2) - round(w_k/4) ]))**2)
		k_e = round(log(R_e,q))
		err = abs(k - (k_s+k_e))
		if err<err_in_solution:
			err_in_solution = err
			w_e = w_k
			Reps_e = R_e
			k_e_good = k_e

	assert(err_in_solution<=1)
	k_e = k_e_good

	k = k_s + k_e

	# upper-lists, we have 4 of them
	L2_error = max(1, multinom(k, [ceil(w_e/8),ceil(w_e/8), k - 2*ceil(w_e/8)]))
	L2 = L2_original*L2_error

	Reps = Reps_s*Reps_e

	L1 = max(1,L2*L2 / Reps)
	L0 = max(1,L1**2 / 2**(floor(n-k)))

	runtime_list = [log(L2_original,2).n(), log(L2_error,2).n(), log(L2,2).n(), log(L1,2).n(), log(L0,2).n()]
	runtime = max(runtime_list)
	lvl = runtime_list.index(runtime)


	return runtime, lvl, runtime_list, k_e, k_s, k

"""
print("--------- Rep-0 with reps for error---------")
for param in all_NTRU_ieee:
	res1 = Rep0Error(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', res1)
"""
# ------REP0-1-----
def Rep1Error(params, eps_s, eps_e):
	n = params['n']
	q = params['q']
	w = params['w']
	we = params['dg']

	Reps_s = (multinom(round(w/2), [round(w/4),round(w/2) - round(w/4) ]))**2 * multinom(n - 2*round(w/2), [eps_s, eps_s, n - 2*round(w/2) - 2*eps_s ] )
	L2_s = multinom(ceil(n/2), [ceil(w/8)+eps_s/2,ceil(w/8)+eps_s/2, ceil(n/2) - 2*ceil(w/8) - eps_s])

	#print(log(Reps_s,2).n())

	k_s = round(log(Reps_s, q))

	#brute-force search for k (# of coordinates we merge on) that satisfies
	# k = k_s + log(Reps_e,2), where Reps_e is the number of representations
	# for the error on k coordinates

	err_in_solution = 100
	w_e = 0
	k_good = 0
	for k in range(k_s, round(n/2), 2):
		# expected error weight on k
		w_k = round(k*we/n)
		if w_k %2 == 1:
			w_k+=1

		R_e = max(1, (multinom(round(w_k/2), [round(w_k/4),round(w_k/2) - round(w_k/4) ]))**2 * multinom(k - 2*round(w_k/2), [eps_e, eps_e, k - 2*round(w_k/2) - 2*eps_e ] ))
		k_e = round(log(R_e,q))
		err = abs(k - (k_s+k_e))
		if err<err_in_solution:
			err_in_solution = err
			w_e = w_k
			Reps_e = R_e
			k_e_good = k_e

	assert(err_in_solution<=1)

	#k_e = k_good
	#print(k_good, err_in_solution, log(Reps_e, 2).n())
	k = k_s + k_e_good


	Cost = {}
	if k/2 < 2*ceil(w_e/8)+eps_e:
		Cost["rt"] = float("infinity")
		return Cost

	# upper-lists, we have 4 of them
	L2_error = max(1, multinom(k, [ceil(w_e/8)+eps_e/2,ceil(w_e/8)+eps_e/2, k - 2*ceil(w_e/8)-eps_e]))

	L2 = L2_s*L2_error

	#print(log(L2_s,2).n(), log(L2_error,2).n(), w_e)

	Reps = Reps_s*Reps_e

	#print('Reps:', log(Reps,2).n(), 'L2:', log(L2,2).n())

	L1 = max(1,L2*L2 / Reps)
	L0 = max(1,L1**2 / 2**(n-k))


	runtime_list = [log(L2_s,2).n(), log(L2_error,2).n(), log(L2,2).n(), log(L1,2).n(), log(L0,2).n()]
	runtime = max(runtime_list)
	lvl = runtime_list.index(runtime)


	Cost["rt"] = runtime
	Cost["k"] = [k, k_s, k_e_good]
	Cost["list"] = runtime_list
	Cost["lvl"] = lvl
	Cost["we"] = w_e


	return Cost

"""
print("--------- Rep-1 of depth 3 with reps for error---------")
for param in all_param_sets:
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', Rep1Error(param, 0, 0))

	eps_s_range = [i for i in range(0,20,2)]
	eps_s_range.reverse()

	eps_e_range = [j for j in range(0,20,2)]
	eps_e_range.reverse()

	min_cost = {}
	min_cost["rt"] = float("infinity")

	for i in eps_s_range:
		for j in eps_e_range:
			cost = Rep1Error(param, i, j)
			if cost["rt"] < min_cost["rt"]:
				min_cost = cost

#print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', min_cost, min_eps)

"""


#brute-force search for k (# of coordinates we merge on) that satisfies
# k = k_s + log(Reps_e,2), where Reps_e is the number of representations
# for the error on k coordinates
def find_k(n, q, w, k_s, eps_e, k_zeroized):
	err_in_solution = 100

	k_e = 0
	Reps_e = 1
	w_e = 0


	for k in range(max(2,k_zeroized+2), k_s*2, 2):
		# expected error weight on k for lvl

		k_remaining = k - k_zeroized
		w_k = round(k_remaining*w/n)
		#if w_k %2 == 1:
		#	w_k+=1

		R_e = max(1, (multinom(round(w_k/2), [round(w_k/4),round(w_k/2) - round(w_k/4) ]))**2 * multinom(k_remaining - 2*round(w_k/2), [eps_e, eps_e, k_remaining - 2*round(w_k/2) - 2*eps_e ] ))
		k_e = round(log(R_e,q))
		err = abs(k - (k_s+k_e))
		if err<err_in_solution:
			err_in_solution = err
			w_e = w_k
			Reps_e = R_e
			k_e_good = k_e
			k_good = k

	return k_e_good, k_good, Reps_e, w_e

def Rep1ErrorDepth4_version2(params, eps_s, eps_e):
	n = params['n']
	q = params['q']
	w = params['w']
	w_e = params['dg']

	n_lvl = len(eps_s)
	ws = vector(RR, n_lvl+1)
	we = vector(RR, n_lvl+1)

	ws[0] = round(w/2)
	we[0] = round(w_e/2)

	for i in range(1,n_lvl+1):
		ws[i] = round(ws[i-1]/2) + eps_s[i-1]
		we[i] = round(we[i-1]/2) + eps_e[i-1]

	#print('ws:', ws)


	Reps = vector(RR, n_lvl)   # number of representations
	k    = vector(ZZ, n_lvl+1) #length of window on which we match
	ke   = vector(ZZ, n_lvl+1)
	ks   = vector(ZZ, n_lvl+1)
	S 	 = vector(RR,n_lvl)    # search space
	L 	 = [0]*(n_lvl+2)#vector(RR,n_lvl+2)  #list sizes
	T 	 = vector(RR,n_lvl+2)  #runtimes
	we   = vector(RR, n_lvl+1) # weights of e on ke

	Cost = {}

	# compute representation and merging windows from bottom to top
	i_reverse_range = [i for i in range(n_lvl)]
	i_reverse_range.reverse()
	for i in i_reverse_range:
		Reps_s = multinom(ws[i], [round(ws[i]/2), ws[i]-round(ws[i]/2)])**2 * multinom(n - 2*ws[i], [eps_s[i], eps_s[i], n - 2*ws[i] - 2*eps_s[i] ] )
		k_s = round(log(Reps_s, q))

		k_e, k_, Reps_e, w_e = find_k(n, q, w_e, k_s, eps_e[i], k[n_lvl-i-1]) #k[0]=0
		#print('Reps_e:', log(Reps_e,2).n(), 'Reps_s:', log(Reps_s,2).n())
		Reps[i] = Reps_s * Reps_e
		we[n_lvl-i] = w_e #weight of the error on the coordinates we enumerate (k_ - k[n_lvl-i-1])
		ke[n_lvl-i] = k_e
		ks[n_lvl-i] = k_s
		k[n_lvl-i] = k_



	# compute the search spaces and the list sizes
	# start with top-most lists (mitm enumeration of s, enumeration of e on k[1])
	Ltop_s = multinom(ceil(n/2), [ceil(ws[n_lvl]/2), ceil(ws[n_lvl]/2), ceil(n/2) - 2*ceil(ws[n_lvl]/2)])
	#k[1] is the length of the first merge window
	if k[1]<=2*(round(we[1]/4)+eps_e[n_lvl-1]):
		#print('k[1] is too small:', k[1])
		Cost["rt"] = float("infinity")
		return Cost
	Ltop_error = max(1, multinom(k[1], [round(we[1]/4)+eps_e[n_lvl-1],round(we[1]/4)+eps_e[n_lvl-1], k[1] - 2*round(we[1]/4)-2*eps_e[n_lvl-1]]))
	L[0] = Ltop_s * Ltop_error

	#print('Ltop_s:', log(Ltop_s,2).n())

	T[0] = round(log(L[0],2))

	# compute the search spaces and the list sizes from top to bottom
	# except the outer ones
	for i in i_reverse_range:
		S[i] = multinom(n, [ws[i+1], ws[i+1], n - 2*ws[i+1]])
		if i>0:
			if k[i+1]-k[i]<=2*(round(we[i+1]/2)+eps_e[i-1]):
				#print('k[i+1]-k[i] is too small:', k[i+1]-k[i], we[i+1], eps_e[i-1])
				Cost["rt"] = float("infinity")
				return Cost
			S[i] = S[i]*multinom(k[i+1]-k[i], [round(we[i+1]/2)+eps_e[i-1], round(we[i+1]/2)+eps_e[i-1],k[i+1]-k[i] -2*round(we[i+1]/2)-2*eps_e[i-1]])
		L[n_lvl-i] = S[i] / Reps[i]
	#L[0] = sqrt(S[n_lvl-1]) #should be ~ L[0] = Ltop_s * Ltop_error

	L[n_lvl+1] = L[n_lvl]**2/(2**(n-k[n_lvl]))

	#print('R:',[log(Reps[i],2).n() for i in range(len(Reps))])
	#print('S:', [log(S[i],2).n() for i in range(len(S))])
	#print('L:', [log(L[i],2).n() for i in range(len(L))])

	for i in range(1, n_lvl+1):
		#max{expected (unfiltered) output-size, input-size)}
		T[i] = round( log( max( ( ( L[i-1]**2 )/(q**(k[i]-k[i-1]))),(L[i-1]) ), 2))

	T[n_lvl+1] = max(1, round(log(L[n_lvl+1],2)))

	runtime = max(T)
	lvl = list(T).index(runtime)

	Cost["rt"] = runtime
	Cost["runtimes"] = T
	Cost["lists"] =[log(L[i],2).n() for i in range(len(L))]
	Cost["k"] = k
	Cost["ke"] = ke
	Cost["ks"] = ks
	Cost["we"] = we
	Cost["lvl"] = lvl

	return Cost


"""
print("--------- Rep-1 of depth 4 with reps for error---------")
for param in all_NTRU_ieee:
	#print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', Rep1ErrorDepth4_version2(param, [0], [0]))

	eps_s_range1 = [i for i in range(0,20,2)]
	eps_s_range1.reverse()

	eps_e_range1 = [j for j in range(0,10,2)]
	eps_e_range1.reverse()

	min_cost = {}
	min_cost["rt"] = float("infinity")

	for i1 in eps_s_range1:
		for j1 in eps_e_range1:
			eps_s_range2 = [i for i in range(0,max(i1,1), 2)]
			eps_s_range2.reverse()

			eps_e_range2 = [j for j in range(0,max(j1,1), 2)]
			eps_e_range2.reverse()

			for i2 in eps_s_range2:
				for j2 in eps_e_range2:
					#cost = Rep1ErrorDepth4(param, [i1,i2], [j1,j2])
					cost = Rep1ErrorDepth4_version2(param, [i1,i2], [j1,j2])
					#print([i1,i2], [j1,j2], cost["rt"])
					if cost["rt"] < min_cost["rt"]:
						#print(min_cost["rt"])
						min_cost = cost
						min_eps_s = [i1,i2]
						min_eps_e = [j1,j2]

	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', min_cost, min_eps_s, min_eps_e)
	#print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', Rep1ErrorDepth4_version2(param, [0,0], [0,0]))

#print(Rep1ErrorDepth4_version2(paramset_NTRU1, [10,8], [2,0]))
#print(Rep1ErrorDepth4_version2(paramset_GLP, [12,4], [0,0]))

"""

print("--------- Rep-1 of depth 5 with reps for error---------")
for param in all_NTRU_ieee:
	eps_s_range1 = [i for i in range(18,36,2)]
	eps_s_range1.reverse()

	eps_e_range1 = [j for j in range(0,12,2)]
	eps_e_range1.reverse()

	min_cost = {}
	min_cost["rt"] = float("infinity")

	for i1 in eps_s_range1:
		for j1 in eps_e_range1:

			eps_s_range2 = [i2 for i2 in range(0,max(1,20), 2)]
			eps_s_range2.reverse()

			eps_e_range2 = [j2 for j2 in range(0,max(1,j1-4), 2)]
			eps_e_range2.reverse()

			for i2 in eps_s_range2:
				for j2 in eps_e_range2:

					eps_s_range3 = [i3 for i3 in range(0,max(1,i2), 2)]
					eps_s_range3.reverse()

					eps_e_range3 = [j3 for j3 in range(0,max(1,j2-2), 2)]
					eps_e_range3.reverse()

					for i3 in eps_s_range3:
						for j3 in eps_e_range3:

							cost = Rep1ErrorDepth4_version2(param, [i1,i2,i3], [j1,j2,j3])
							if cost["rt"] < min_cost["rt"]:
								print('new min:', cost["rt"])
								min_cost = cost
								min_eps_s = [i1,i2,i3]
								min_eps_e = [j1,j2,j3]
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', min_cost, min_eps_s, min_eps_e)

"""
print("--------- Rep-1 of depth 6 with reps for error---------")
for param in all_param_sets:
	eps_s_range1 = [i for i in range(0,20,2)]
	eps_s_range1.reverse()

	eps_e_range1 = [j for j in range(0,12,2)]
	eps_e_range1.reverse()

	min_cost = {}
	min_cost["rt"] = float("infinity")

	for i1 in eps_s_range1:
		for j1 in eps_e_range1:

			eps_s_range2 = [i2 for i2 in range(0,max(1,i1), 2)]
			eps_s_range2.reverse()

			eps_e_range2 = [j2 for j2 in range(0,max(1,j1), 2)]
			eps_e_range2.reverse()

			for i2 in eps_s_range2:
				for j2 in eps_e_range2:

					eps_s_range3 = [i3 for i3 in range(0,max(1,i2), 2)]
					eps_s_range3.reverse()

					eps_e_range3 = [j3 for j3 in range(0,max(1,j2), 2)]
					eps_e_range3.reverse()

					for i3 in eps_s_range3:
						for j3 in eps_e_range3:

							eps_s_range4 = [i4 for i4 in range(0,max(1,i3), 2)]
							eps_s_range4.reverse()

							eps_e_range4 = [j4 for j4 in range(0,max(1,i3), 2)]
							eps_e_range4.reverse()

							for i4 in eps_s_range4:
								for j4 in eps_e_range4:
									cost = Rep1ErrorDepth4_version2(param, [i1,i2,i3,i4], [j1,j2,j3,j4])
									if cost["rt"] < min_cost["rt"]:
										print('new min:', cost["rt"],[i1,i2,i3,i4], [j1,j2,j3,j4])
										min_cost = cost
										min_eps_s = [i1,i2,i3,i4]
										min_eps_e = [j1,j2,j3,j4]
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', min_cost, min_eps_s, min_eps_e)
"""
