

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
		res*=binomial(n_, c[i])
		n_ = n_ - c[i]
	return res

# ------REP0------
def Rep0(params):
	n = params['n']
	q = params['q']
	w = params['w']

	min_runtime = float("infinity")

	# upper-lists, we have 4 of them
	L2 = multinom(ceil(n/2), [ceil(w/8),ceil(w/8), ceil(n/2) - 2*ceil(w/8)])
	T2 = L2

	# number of representations
	Reps1 = (multinom(round(w/2), [round(w/4),round(w/2) - round(w/4) ]))**2

	k = floor(log(Reps1,2) / (log(q,2) - 0.5*log(3,2)))
	if k % 2 == 1:
		k = k-1

	# search for optimal B for NN search
	for B in range(3, 50):

		#nReps = round((1 - 1/B)^(-k/2)) # (1 - 1/3)^(-k/2)
		nReps = ((B/(B-1))**k)*n
		#L1 = max(1,L2*L2 / Reps1)
		#L1 = multinom(n, [ceil(w/4), ceil(w/4), n - 2*ceil(w/4)]) / Reps1

		#expected output in one run
		Bucket_size = L2 *(1/q)^(k/2)*(B/q)^(k/2)
		T1 = q^(k/2) * (ceil(q/B))^(k/2)*Bucket_size**2 #number of labels X naive search inside each bucket

		L1 = L2^2 * 3^(k/2) / q^k
		T0 = max(1, L1**2 / 2**(floor(n-k)))
		L0 = L1**2 / 2**(floor(n-k))

		runtime_list = [T2, T1*nReps, nReps*T0]
		runtime = max(runtime_list)
		lvl = runtime_list.index(runtime)

		if runtime<min_runtime:
			min_runtime = runtime
			min_L = [log(L2,2).n(), log(L1,2).n(), log(L0,2).n()]
			min_T = [log(T2,2).n(), log(T1*nReps,2).n(), log(T0,2).n()]
			min_nReps = nReps
			min_k = k
			min_lvl = lvl
			min_B = B

	Cost = {}
	Cost["runtime"] = log(min_runtime,2).n()
	Cost["T"] = min_T
	Cost["list"] = min_L
	Cost["NN_nReps"] = log(min_nReps,2).n()
	Cost["k"] = min_k
	Cost["lvl"] = min_lvl
	Cost["B"] = min_B


	return Cost
"""
print("--------- Rep-0 ---------")
for param in all_NTRU_ieee:
	runtime = Rep0(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime)
"""


# ------REP1 (same depth as Rep0)------
def Rep1_depth3(params, eps):
	n = params['n']
	q = params['q']
	w = params['w']

	ws = round(w/2) + eps

	# upper-lists, we have 4 of them
	L2 = multinom(ceil(n/2), [ceil(w/8)+eps/2,ceil(w/8)+eps/2, ceil(n/2) - 2*ceil(w/8) - eps])
	T2 = L2

	Reps10 = multinom(round(w/2), [round(w/4), round(w/2)-round(w/4)])**2 * multinom(n - w, [eps, eps, n - w - 2*eps ] )
	k = round(log(Reps10,2) / (log(q,2) - 0.5*log(3,2)))
	if not k%2==0:
		k = k-1
	#print(log(Reps10,2).n(), k)


	L1 = L2^2 * 3^(k/2) / q^k
	L0 = L1**2 / 2**(floor(n-k))

	#print([log(L0,2).n(), log(L1,2).n(), log(L2,2).n()])

	min_runtime = float("infinity")
	for B in range(3, 30):


		#rho = log(1 - 1/B) / log(B/q)
		#nReps = ceil(L2^rho)
		nReps = ((B/(B-1))**k)*n

		Bucket_size = max(1,ceil(L2 *(1/q)^(k/2)*((B/q))^(k/2)))

		#print(B, rho.n(), nReps, log(Bucket_size,2).n())
		T1 = q^(k/2) * ((q/B))^(k/2)*Bucket_size**2

		#L1 = max(1,L2*L2 / Reps10)

		T0 = max(1, L1**2 / 2**(floor(n-k)))

		runtime_list = [T2, T1*nReps, nReps*T0]
		runtime = max(runtime_list)
		lvl = runtime_list.index(runtime)

		if runtime<min_runtime:
			min_runtime = runtime
			min_L = [log(L0,2).n(), log(L1,2).n(), log(L2,2).n()]
			min_T = [log(nReps*T0,2).n(), log(T1*nReps,2).n(), log(T2,2).n()]
			min_bucketsize = log(Bucket_size,2).n()
			min_nReps = nReps
			min_k = k
			min_lvl = lvl
			min_B = B


	Cost = {}
	Cost["runtime"] = round(log(min_runtime,2))
	Cost["T"] = min_T
	Cost["list"] = min_L
	Cost["NN_nReps"] = log(min_nReps,2).n()
	Cost["k"] = min_k
	Cost["lvl"] = min_lvl
	Cost["B"] = min_B

	return Cost



def R1_depth3_eps(param):
	eps1_range = [i for i in range(0,30,2)]
	eps1_range.reverse()

	min_cost = {}
	min_cost["runtime"] = float("infinity")

	for i in eps1_range:
		cost = Rep1_depth3(param, i)
		if cost["runtime"] < min_cost["runtime"]:
			min_cost = cost
			min_eps = i
	return min_cost, min_eps



def R1(params, eps):
	n = params['n']
	q = params['q']
	w = params['w']


	n_lvl = len(eps)

	ws = vector(RR, n_lvl+1)
	ws[0] = round(w/2)
	for i in range(1,n_lvl+1):
		ws[i] = round(ws[i-1]/2) + eps[i-1]

	#print('w:', ws)

	Reps = vector(RR, n_lvl) #[0]*(n_lvl)
	r    = vector(RR, n_lvl+1)#[0]*(n_lvl+1)#
	S 	 = vector(RR,n_lvl)#[0]*(n_lvl)#


	L 	= [0]*(n_lvl+2) #list sizes
	L_  = [0]*(n_lvl+2)
	T 	= [0]*(n_lvl+2) #runtimes
	k 	= [0]*(n_lvl+1) #length of windows we merge on
	Bucket_size = [0]*n_lvl
	nRepeat = [0]*n_lvl

	Cost = {}
	for i in range(n_lvl):
		#number of representations
		assert n-2*ws[i] - 2*eps[i] > 1, "too large eps and ws"
		Reps[i] = multinom(ws[i], [round(ws[i]/2), ws[i]-round(ws[i]/2)])**2 * multinom(n - 2*ws[i], [eps[i], eps[i], n - 2*ws[i] - 2*eps[i] ] )
		#print(log(Reps[i],2).n())
		#search spaces
		S[i] 	= multinom(n, [ws[i+1], ws[i+1], n - 2*ws[i+1]])
		#k[i] = round(2^(i+1)*log(Reps[i],2)/( (2^(i+1)-1)*log(q,2)-log(3,2)))
		#L[i+1] = ceil(S[i] / Reps[i])
		k[i] = ceil(log(Reps[i],2)/(log(q,2) - (0.5)^(i+1)*log(3,2)) )
		#if((k[i]+1) % 2^(i+1)==0):
		#	k[i] = k[i]+1
		while not (k[i] % 2^(i+1) == 0):
			k[i]-= 1

		if(k[i]<=0):
			Cost["runtime"] = float("infinity")
			return Cost
	#print(log(Reps[0],2).n(), k[0])
	#L[n_lvl+1] = sqrt(S[n_lvl - 1]) # size of the top-most lists, close to the value below
	L[n_lvl+1] = multinom(ceil(n/2), [ceil(ws[n_lvl]/2), ceil(ws[n_lvl]/2), ceil(n/2) - 2*ceil(ws[n_lvl]/2)])

	#k1 = log(Reps[1],2)/(log(q,2) - (0.5)^(2)*log(3,2))
	#print(k1.n(), log(L[2],2).n(), log(L[3]^2*(1/q)^(3*k1/4)*(3/q)^(k1/4),2).n(), log(L[3]^2*(1/q)^(3*k[1]/4)*(3/q)^(k[1]/4),2).n() )

	#k0 = log(Reps[0],2)/(log(q,2) - (0.5)*log(3,2))
	#print(k0.n(), log(L[1],2).n(), log(L[2]^2*(1/q)^(k0/2)*(3/q)^(k0/2),2).n(),log(L[2]^2*(1/q)^((k[0]-k[1])/2)*(3/q)^((k[0]-k[1])/2),2).n())
	#print([log(L[i],2).n() for i in range(len(L))])

	min_runtime = float("infinity")

	T[n_lvl+1] = L[n_lvl+1]



	if n_lvl<=2:
		C = 35
	else:
		C = 12
	maxrange = C^n_lvl
	B = [0]*n_lvl
	for B_ in range(3, maxrange):
		tmp = B_
		for i in range(n_lvl):
			if i == n_lvl - 1:
				B[i] = max(3,tmp)
			else:
				B[i] = max(3, tmp // C^(n_lvl-i-1))
				tmp = tmp % C^(n_lvl-i-1)
			assert(B[i]<C)

		# time to construct middle lists
		for i in range(n_lvl):
			ind = n_lvl - i - 1
			#rho = log(1 - 1/B[ind]) / log(B[ind]/q)
			#nRepeat[ind] = ceil(L[n_lvl-i+1]^rho)

			nRepeat[ind] = (B[ind]/(B[ind]-1))**(k[ind]-k[ind+1])*n
			pow_q  = (k[ind]-k[ind+1])*(2^(ind+1)-1) / (2^(ind+1))
			pow_bq = (k[ind]-k[ind+1])/ (2^(ind+1))

			L[n_lvl-i] = ceil(L[n_lvl-i+1]^2*(1./q)^(pow_q)*(3./q)^pow_bq)

			Bucket_size[ind] = max(1,ceil(L[n_lvl+1-i] *(1/q)^(pow_q)*(B[ind]/q)^pow_bq))


			T[ind+1] = (nRepeat[ind]*q^pow_q * (q/B[ind])^pow_bq) * Bucket_size[ind]**2
			#print(B, ind+1, log(T[ind+1],2).n(), nRepeat[i])

		#print(B, rho.n(), log(T[2],2).n(), log(T[1],2).n())

		L[0] = L[1]**2/(2**(n-k[0]))
		T[0] = max(L[1], L[0])


		runtime = max(T)
		lvl = T.index(runtime)
		if runtime<min_runtime:
			min_runtime = runtime
			min_runtime_list = copy(T)
			min_list = copy(L)
			min_bucketsize = copy(Bucket_size)
			min_nReps = copy(nRepeat)
			min_k = copy(k)
			min_lvl = lvl
			min_B = copy(B)
			#print(min_B, round(log(min_runtime,2)), eps)


	Cost = {}
	Cost["runtime"] = (log(min_runtime,2)).n()
	Cost["T"] = [log(min_runtime_list[i],2).n() for i in range(len(T))]
	Cost["L"] = [log(min_list[i],2).n() for i in range(len(L))]
	Cost["Repeat"] = [log(min_nReps[i],2).n() for i in range(len(min_nReps))]
	Cost["bucket sizes"] = [log(min_bucketsize[i],2).n() for i in range(len(min_bucketsize))]
	Cost["R"] = [log(Reps[i],2).n() for i in range(len(Reps))]
	Cost["lvl"] = min_lvl
	Cost["k"] = min_k
	Cost["B"] = min_B

	return Cost

def R1_depth4_eps(param):
	eps1_range = [i for i in range(6)]
	eps1_range.reverse()
	min_runtime = {}
	min_runtime["runtime"] = float("infinity")
	eps = []
	for i in eps1_range:
		eps2_range = [j for j in range(4)]
		eps2_range.reverse()
		for j in eps2_range:
			cost = R1(param, [i,j])
			#print([i,j], cost["T"])
			if cost["runtime"]<min_runtime["runtime"]:
				#print('new runtime:', cost["runtime"], cost["T"], [i,j])
				min_runtime = cost
				eps = [i,j]
	return min_runtime, eps

"""
print("--------- Rep-1, depth 4 ---------")
for param in all_param_sets:
	runtime, eps = R1_depth4_eps(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime, eps)
"""
"""
print("--------- Rep-1, depth 4 ---------")
for param in all_NTRU_ieee:
	runtime, eps = R1_depth4_eps(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime, eps)
"""
# for tests
def R1_depth3_eps_2(param):
	eps1_range = [i for i in range(0,22,2)]
	eps1_range.reverse()
	min_runtime = {}
	min_runtime["runtime"] = float("infinity")
	eps = []
	for i in eps1_range:
		cost = R1(param, [i])
		if cost["runtime"]<min_runtime["runtime"]:
			#print('new runtime:', cost["runtime"], cost["T"], [i,j])
			min_runtime = cost
			eps = [i]
	return min_runtime, eps
"""
print("--------- Rep-1, depth 3 ---------")
for param in all_param_sets:
	runtime, eps = R1_depth3_eps(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime, eps)
"""

def R1_depth5_eps(param):
	eps1_range = [i for i in range(25,40)]
	eps1_range.reverse()
	min_runtime = {}
	min_runtime["runtime"] = float("infinity")
	for i in eps1_range:
		eps2_range = [j for j in range(20)]
		eps2_range.reverse()
		for j in eps2_range:
			eps3_range = [k for k in range(9)]
			eps3_range.reverse()
			for k in eps3_range:
				cost = R1(param, [i,j,k])
				if cost["runtime"]<min_runtime["runtime"]:
					print('new runtime:', cost["runtime"], cost["T"], [i,j,k])
					min_runtime = cost
					eps = [i,j,k]
	return min_runtime, eps

"""
print("--------- Rep-1, depth 5 ---------")
for param in all_NTRU_ieee:
	runtime, eps = R1_depth5_eps(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime, eps)
"""

def cold_boot(param, alg = R1_depth3_eps):
	n = param['n']

	ws = ceil((1.0+0.1)/(2*100)*2*n)
	if ws%2 == 1:
		ws+=1
	param['w'] = ws
	return alg(param), ws

def cold_boot_guessing(param, alg=R1_depth3_eps):
	n = param['n']

	ws = ceil((1.0+0.1)/(2*100)*2*n)
	if ws%2 == 1:
		ws+=1
	param['w'] = ws

	min_overall = {}
	min_overall["runtime_overall"] = float("infinity")

	for c in range(10, 100, 5):
		p0 = log(binomial(n-c, ws)/binomial(n, ws),2).n() #probability that on c coordinates the weight is 0
		param['n'] = n-c
		runtime1, eps = alg(param)
		runtime_overall = runtime1["runtime"]-p0
		if runtime_overall<min_overall["runtime_overall"]:
			min_overall = runtime1
			min_overall["c"] = c
			min_overall["runtime_overall"] =  runtime_overall
			min_overall["eps"] = eps
	return min_overall



for param in all_param_sets:
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', cold_boot(param, alg=R1_depth4_eps))
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', cold_boot_guessing(param, alg=R1_depth4_eps))

def R1_depth6_eps(param):
	eps1_range = [i for i in range(40,45)]
	eps1_range.reverse()
	min_runtime = {}
	min_runtime["runtime"] = float("infinity")
	for i in eps1_range:
		eps2_range = [j for j in range(25)]
		eps2_range.reverse()
		for j in eps2_range:
			eps3_range = [k for k in range(15)]
			eps3_range.reverse()
			for k in eps3_range:
				eps4_range = [l for l in range(15)]
				eps4_range.reverse()
				for l in eps4_range:
					cost = R1(param, [i,j,k,l])
					if cost["runtime"]<min_runtime["runtime"]:
						print('new runtime:', cost["runtime"], cost["T"], [i,j,k,l])
						min_runtime = cost
						eps = [i,j,k,l]
		print('i:', i)
	return min_runtime, eps
"""
print("--------- Rep-1, depth 6 ---------")
for param in all_param_sets[:1]:
	runtime, eps = R1_depth6_eps(param)
	print('[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], ' | ', runtime, eps)
"""
