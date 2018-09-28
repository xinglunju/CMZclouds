
import time
import numpy as np

start = time.clock()

def clustermass(N=4, dr=1.0):
	'''
	Input:
	N -- number of high-mass protostars, from 1 to N, should be an integer
	dr -- H2O maser (or any other SF tracer) detection rate, between 0 and 1
	# Some references: Szymczak 2005 -- 0.52
	# Cyganowski 2013 -- 0.68
	# Titmarsh 2016 -- 0.50
	'''
	# Define the IMF, following Kroupa 2001, MNRAS, 322, 231
	mtop = 150.
	m  = np.arange(0.01,mtop,0.01)
	lowm = np.arange(0.01,0.08,0.01)[0:7]
	intm = np.arange(0.08,0.50,0.01)
	higm = np.arange(0.50,mtop,0.01)
	# Generate the IMF based on stellar masses
	lowm_pdf = lowm**-0.3 * 12.5
	intm_pdf = intm**(-1.3)
	higm_pdf = higm**(-2.3) * 0.5
	pm = np.concatenate((lowm_pdf,intm_pdf,higm_pdf))
	# Normalization
	pm /= pm.sum()
	meanM = np.sum(m*pm)

	# This is the prior, choose the range of cluster masses
	# based on the number of high-mass protostars
	#step = 10.
	#Mcluster = np.arange(1e1,2e3,step)
	ends = np.logspace(np.log10(10.),np.log10(2000.),300)
	Mcluster = (ends[1:] + ends[:-1])/2.
	dm = np.diff(ends)

	dat = file('HMdet_'+str(N)+'_dr_'+str(dr)+'_log.dat', 'ab')

	# Likelihood for each Mcluster, and high-mass protostars from 1 to 10
	#likelihood = np.zeros((len(Mcluster),N))
	for ind in np.arange(0,len(Mcluster)):
		mc = Mcluster[ind]
		n_runs = 10000
		match = np.zeros(N)
		tolerance = dm[ind]/2.
		for k in np.arange(n_runs):
			m_stars = np.array([])
			mtot = 0
			while (mtot - mc) < tolerance:
				morestars = max(int(np.ceil((mc+tolerance-mtot) / meanM)), 1)
				m_morestars = np.random.choice(m, p=pm, size=morestars)
				m_stars = np.concatenate([m_stars, m_morestars])
				mtot = m_stars.sum()
			if (mtot - mc) > tolerance:
				mcum = m_stars.cumsum()
				last_ind = np.argmin(np.abs(mcum - mc)) + 1
				m_stars = m_stars[:last_ind]
				mtot = m_stars.sum()
			if (mtot - mc) < -tolerance:
				k -= 1
				continue
			n_high = np.count_nonzero(m_stars >= 8)
			if dr < 1:
				cc = np.random.choice([1,0], p=[dr,1-dr], size=n_high)
				n_high = cc.sum()
			if 1 <= n_high <= N:
				match[n_high-1] += 1
		#likelihood[ind, :] = match / n_runs
		#np.savetxt(dat, np.c_[Mcluster,likelihood], delimiter='\t')
		np.savetxt(dat, np.concatenate(((mc,), match)).reshape(1,N+1), delimiter='\t')
		print "Finished running Mcluster of", mc
	
	dat.close()

clustermass(N=10,dr=1.0)

elapsed = (time.clock() - start)
print 'Time used: {0:0.0f} seconds, or {1:0.1f} minutes.'.format(elapsed, elapsed/60.)
