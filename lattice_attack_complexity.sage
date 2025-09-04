

load("estimator.py")

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

def lattice_optimistic(params):
	n = params['n']
	q = params['q']
	w = params['w']
	dg = params['dg']

	sd = sqrt(dg / n)
	alpha = sqrt(2*pi)*sd/RR(q)
	m = n
	secret_distribution = ((-1, 1), w) #w compinents are (+/-1), the rest is 0
	success_probability = 0.99

	#optimistic quantum sieve:
	#reduction_cost_model =  lambda beta, d, B: ZZ(2)**RR(0.265*beta)

	#quantum enumeration:
	#reduction_cost_model =  lambda beta, d, B: ZZ(2)**RR(0.5*(0.187*beta*log(beta,2)-1.019*beta + 16.1))

	#classical optimistic:
	reduction_cost_model =  lambda beta, d, B: ZZ(2)**RR(0.292*beta+16.4)
	primald = partial(drop_and_solve, primal_usvp, postprocess=False, decision=False)
	return primald(n, alpha, q, secret_distribution=secret_distribution, m=m,  success_probability=success_probability, reduction_cost_model=reduction_cost_model)


#print('---------------------- Lattices  ---------------------- ')
#for param in all_param_sets:
#	lattice_complexity = lattice_optimistic(param)
#	print('classical lattice:[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], lattice_complexity)


print('---------------------- Lattices  ---------------------- ')
for param in all_param_sets:
	lattice_complexity = lattice_optimistic(param)
	print('classical lattice:[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], lattice_complexity)

# print('---------------------- Lattices cold boot  ---------------------- ')
# for param in all_param_sets:
# 	n = param['n']

# 	ws = ceil((1.0+0.1)/(2*100)*2*n)
# 	if ws%2 == 0:
# 		ws+=1
# 	param['w'] = ws

# 	lattice_complexity = lattice_optimistic(param)
# 	print('classical lattice:[n =', param['n'], ' q = ', param['q'], 'w = ', param['w'], lattice_complexity)
