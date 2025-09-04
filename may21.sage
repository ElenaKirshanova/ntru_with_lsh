from sage.all import *


sparceLWE1 = {'n': 512, 'q': 2**12, 'w': 11, 'eta': 2}  #stddev = sqrt(eta/2.0)
sparceLWE2 = {'n': 1024, 'q': 2**26, 'w': 11, 'eta': 2} #stddev = sqrt(eta/2.0)
all_sparse = [sparceLWE1, sparceLWE2]

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
				 paramset_BLISS,paramset_GLP, 
				 NTRU_ieee401, NTRU_ieee449, NTRU_ieee677, NTRU_ieee1087,
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

def H(n,c):
	"""
	Entrppy function
	"""
	assert sum(c) == n, 'bad input to H!'
	res = 0
	for i in range(len(c)):
		res-=(c[i]/n)*log(c[i]/n, 2)

	return res.n()

def binomial_entropy(eta):
	h = 0
	for i in range(-eta, eta+1, 1):
		tmp = binomial(2*eta, eta+i)
		h -= tmp*log(tmp,2)
	return h

def ternary_entropy(w, n):
	h = - w/n *log(w/(2*n),2) - (n-w)/w * log((n-w)/w ,2)
	return h



def Odlyzko(n, w):
	"""
	Odlyzko’s Meet-in-the-Middle (classical) algorithm
	"""
	L = multinom(ceil(n/2), [ceil(w/4),ceil(w/4), ceil(n/2) - 2*ceil(w/4)])
	return ceil(log(L, 2).n())

	
def RepO(params, error_entropy):
	n = params['n']
	q = params['q']
	w = params['w']

	L2 = multinom(ceil(n/2), [ceil(w/8),ceil(w/8), ceil(n/2) - 2*ceil(w/8)])
	Reps1 = (multinom(round(w/2), [round(w/4),round(w/2) - round(w/4) ]))**2

	#intermediate lists, we have 2 of them
	#L1 = multinom(n, [ceil(w/4), ceil(w/4), n - 2*ceil(w/4)]) / Reps1
	L1 = max(1,ceil(L2**2 / Reps1))

	# expected size of the output list
	L0 = max(1,L1**2 / 2**(floor(n-log(Reps1,q))))

	#assert(Reps1>=q)
	
	mem = max(L0, L1, L2)
	guessing = max(1,2**(error_entropy*ceil(log(Reps1, q))))

	runtime = guessing * mem
	runtime_log = log(runtime,2)

	return ceil(log(L2,2)), ceil(runtime_log), ceil(log(mem,2)), ceil(log(guessing,2))


def Rep1(params, eps, error_entropy,  verbose = False):
	"""
	Represent "0 = 1 - 1 / -1+1" (classical)
	"""
	n = params['n']
	q = params['q']
	w = params['w']

	n_lvl = len(eps)
	ws = vector(RR, n_lvl+1)
	ws[0] = round(w/2)

	for i in range(1,n_lvl+1):
		ws[i] = round(ws[i-1]/2) + eps[i-1]


	Reps = vector(RR, n_lvl) #[0]*(n_lvl)
	r    = vector(RR, n_lvl+1)#[0]*(n_lvl+1)#
	S 	 = vector(RR,n_lvl)#[0]*(n_lvl)#
	L 	 = vector(RR,n_lvl+1) #[0]*(n_lvl+1)
	T 	 = vector(RR,n_lvl+2)#[0]*(n_lvl+2)#vector(RR,n_lvl+2)

	for i in range(n_lvl):
		#number of representations
		Reps[i] = multinom(ws[i], [round(ws[i]/2), ws[i]-round(ws[i]/2)])**2 * multinom(n - 2*ws[i], [eps[i], eps[i], n - 2*ws[i] - 2*eps[i] ] )
		#search spaces
		S[i] 	= multinom(n, [ws[i+1], ws[i+1], n - 2*ws[i+1]])
		#list sizes
		#L[i] = log(S[i],2) - log(Reps[i],2)
		L[i] 	= ceil(S[i] / Reps[i])
		r[i] 	= max(1,round(log(Reps[i], q)))


	L[n_lvl] = ceil(sqrt((S[n_lvl - 1]))) #top-level list-size, we have 2^(n_lvl+2) of them
	#L[n_lvl] = 0.5*log(S[n_lvl - 1],2) # on the log-scale
	r[n_lvl] = 0


	#time to create the bottom-most list max(input, expected output)
	T[0] = max(ceil(L[0]), ceil( (L[0]**2)/(2**(n-r[0])) ) )
	#T[0] = max(gamma*L[0], 2*gamma*L[0] - (n-r[0]) )
	T[0] = round(log(T[0],2))

	#time to create the top-most list
	T[n_lvl+1] = L[n_lvl]
	T[n_lvl+1] = round(log(T[n_lvl+1],2))

	for i in range(1, n_lvl+1):
		#max{expected (unfiltered for the weights of s) output-size, input-size)
		T[i] = max( ( ( L[i]**2 )/(q**(r[i-1]-r[i]))),(L[i]) )
		T[i] = round(log(T[i],2))
		#T[i] = max(2*gamma*L[i] - (r[i-1]-r[i])*log(q,2), gamma*L[i])

	mem = ceil(log(max(L),2))
	runtime = max(T)
	#guessingT = round(r[0]*log(3,2))

	guessingT = max(0,error_entropy*r[0])

	maxlvl = [list(T).index(runtime)]

	if verbose:
		print('T:', [log(T[i], 2) for i in range(len(T))])
		print('L:', [log(L[i], 2) for i in range(len(L))])
		print('guessing:', log(guessingT, 2))

	Cost = {}
	#Cost["classical"] = True
	Cost["runtime"] = round(runtime)
	Cost["guessingT"] = guessingT
	Cost["overall"] = round(runtime+guessingT)
	Cost["maxlvl"] = maxlvl

	return Cost

#	search for optimal eps array (2 dim)
def optimal_R1_depth2(param, entropy, alg=Rep1):
	"""
		search for number of additional "1's"
	"""
	eps1_range = [i for i in range(55)]
	eps1_range.reverse()
	min_runtime = {}
	min_runtime["overall"] = float("infinity")
	eps = []
	for i in eps1_range:
		#return ceil(log(L[n_lvl],2)), ceil(log(runtime,2)), guessingT, ceil(log(update,2)), maxlvl
		#return min_N,minRT,guessing,min_gamma,_, maxlvl
		#cost = alg(param, [i], gammainp=gammainp)
		cost = alg(param, [i], entropy)
		if cost["overall"]<min_runtime["overall"]:
			#print('new min:', memory,time, guessing, overtall_time,gamma)
			min_runtime = cost
			eps = [i]
	return min_runtime, eps

for param in all_param_sets:
	n = param['n']
	q = param['q']
	w = param['w']
	dg = param['dg']
	#eta = param['eta']
	runtimeOld = Odlyzko(n, w)

	#entropy = binomial_entropy(eta)
	entropy = ternary_entropy(dg, n)

	#print("entropy:", entropy.n())
	runtimeRep0 = RepO(param, entropy)
	print('Odlyzko:', param, runtimeOld) # ' with rough poly-factors :', ceil(runtimeOld+log(n/2,2)+log(q,2)))
	print('Rep0:', param, runtimeRep0[1], runtimeRep0[1]+log(n,2), runtimeRep0[2])

	runtumeRep1 = optimal_R1_depth2(param, entropy)
	print('Rep1:', param, runtumeRep1)
	print("------------------------------------")