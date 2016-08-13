from ROOT import *
import numpy as np
import PyROOTPlots as PyRPl
import Poisson_file_handler as Poisson_fh
import matplotlib.pylab as plt

import emcee
import triangle
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

def fit_surface_relation_FID808():

	"""Plot EI vs EH and fit
	for FID808, the bolo with surface calibration
    
	Detail:

	Args:

	Returns:
		void

	Raises:
		void
	"""


	# FWIA	FWIB	FWIC	FWID	FWC1	FWC2	VFID	VET	
	# 1.09	0.74	1.01	0.89	1.12	1.80	8	5.50

	#Open event files
	data_types = {"names": ("RUN", "SN", "EC1", "EC2", "EIA", "EIB", "EIC", "EID"), "formats": ("i", "i", "f", "f", "f", "f",  "f", "f")}

	arr_S1Beta = np.loadtxt("../Text_files/S1Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Beta = np.loadtxt("../Text_files/S2Beta_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S1Pb = np.loadtxt("../Text_files/S1Pb_heatremoved.txt", delimiter=",",  dtype=data_types)
	arr_S2Pb = np.loadtxt("../Text_files/S2Pb_heatremoved.txt", delimiter=",",  dtype=data_types)

	arr_EC_S1Beta, arr_EI_S1Beta = arr_S1Beta["EC1"], arr_S1Beta["EIB"]
	arr_EC_S2Beta, arr_EI_S2Beta = arr_S2Beta["EC1"], arr_S2Beta["EID"]
	arr_EC_S1Pb, arr_EI_S1Pb = arr_S1Pb["EC1"], arr_S1Pb["EIB"]
	arr_EC_S2Pb, arr_EI_S2Pb = arr_S2Pb["EC1"], arr_S2Pb["EID"]

	lS1Beta, lS2Beta, lS1Pb, lS2Pb = np.where(arr_EC_S1Beta<40), np.where(arr_EC_S2Beta<40), np.where(arr_EC_S1Pb<40), np.where(arr_EC_S2Pb<40)

	arr_EI_S1Beta, arr_EC_S1Beta = arr_EI_S1Beta[lS1Beta], arr_EC_S1Beta[lS1Beta]
	arr_EI_S2Beta, arr_EC_S2Beta = arr_EI_S2Beta[lS2Beta], arr_EC_S2Beta[lS2Beta]
	arr_EI_S1Pb, arr_EC_S1Pb     = arr_EI_S1Pb[lS1Pb], arr_EC_S1Pb[lS1Pb]
	arr_EI_S2Pb, arr_EC_S2Pb     = arr_EI_S2Pb[lS2Pb], arr_EC_S2Pb[lS2Pb]

	arr_EI_S1Beta, arr_EC_S1Beta = np.array(arr_EI_S1Beta).astype(float), np.array(arr_EC_S1Beta).astype(float)
	arr_EI_S2Beta, arr_EC_S2Beta = np.array(arr_EI_S2Beta).astype(float), np.array(arr_EC_S2Beta).astype(float)
	arr_EI_S1Pb, arr_EC_S1Pb = np.array(arr_EI_S1Pb).astype(float), np.array(arr_EC_S1Pb).astype(float)
	arr_EI_S2Pb, arr_EC_S2Pb = np.array(arr_EI_S2Pb).astype(float), np.array(arr_EC_S2Pb).astype(float)

	gr_S1Beta, gr_S2Beta   = TGraph(len(arr_EI_S1Beta), arr_EC_S1Beta, arr_EI_S1Beta),  TGraph(len(arr_EI_S2Beta), arr_EC_S2Beta, arr_EI_S2Beta)
	gr_S1Pb, gr_S2Pb       = TGraph(len(arr_EI_S1Pb), arr_EC_S1Pb, arr_EI_S1Pb),  TGraph(len(arr_EI_S2Pb), arr_EC_S2Pb, arr_EI_S2Pb)
	gr_QS1Beta, gr_QS2Beta = TGraph(len(arr_Q_S1Beta), arr_EC_S1Beta, arr_Q_S1Beta),  TGraph(len(arr_Q_S2Beta), arr_EC_S2Beta, arr_Q_S2Beta)
	gr_QS1Pb, gr_QS2Pb     = TGraph(len(arr_Q_S1Pb), arr_EC_S1Pb, arr_Q_S1Pb),  TGraph(len(arr_Q_S2Pb), arr_EC_S2Pb, arr_Q_S2Pb)

	PyRPl.process_TGraph(gr_S1Beta, X_title = "Heat", Y_title = "Ion", color=kOrange-3), PyRPl.process_TGraph(gr_S2Beta, X_title = "Heat", Y_title = "Ion", color=kBlue)
	PyRPl.process_TGraph(gr_S1Pb, X_title = "Heat", Y_title = "Ion", color=kBlack), PyRPl.process_TGraph(gr_S2Pb, X_title = "Heat", Y_title = "Ion", color=kGreen+2)

	PyRPl.process_TGraph(gr_QS1Beta, X_title = "Heat", Y_title = "Q", color=kOrange-3), PyRPl.process_TGraph(gr_QS2Beta, X_title = "Heat", Y_title = "Q", color=kBlue)
	PyRPl.process_TGraph(gr_QS1Pb, X_title = "Heat", Y_title = "Q", color=kBlack), PyRPl.process_TGraph(gr_QS2Pb, X_title = "Heat", Y_title = "Q", color=kGreen+2)


	list_graph = [gr_S1Beta, gr_S2Beta, gr_S1Pb, gr_S2Pb]

	cc = TCanvas("cc", "cc")
	cc.Divide(2,2)
	for i in range(4):
		cc.cd(i+1)
		list_graph[i].Draw("A*")
		list_graph[i].Fit("pol4")









# Reproducible results!
np.random.seed(123)

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)

# Plot the dataset and the true model.
xl = np.array([0, 10])
pl.errorbar(x, y, yerr=yerr, fmt=".k")
pl.plot(xl, m_true*xl+b_true, "k", lw=3, alpha=0.6)
pl.ylim(-9, 9)
pl.xlabel("$x$")
pl.ylabel("$y$")
pl.tight_layout()
pl.savefig("line-data.png")

# Do the least-squares fit and compute the uncertainties.
A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
print("""Least-squares results:
    m = {0} ± {1} (truth: {2})
    b = {3} ± {4} (truth: {5})
""".format(m_ls, np.sqrt(cov[1, 1]), m_true, b_ls, np.sqrt(cov[0, 0]), b_true))

# Plot the least-squares result.
pl.plot(xl, m_ls*xl+b_ls, "--k")
pl.savefig("line-least-squares.png")

# Define the probability function as likelihood * prior.
def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

# Find the maximum likelihood value.
chi2 = lambda *args: -2 * lnlike(*args)
result = op.minimize(chi2, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]
print("""Maximum likelihood result:
    m = {0} (truth: {1})
    b = {2} (truth: {3})
    f = {4} (truth: {5})
""".format(m_ml, m_true, b_ml, b_true, np.exp(lnf_ml), f_true))

# Plot the maximum likelihood result.
pl.plot(xl, m_ml*xl+b_ml, "k", lw=2)
pl.savefig("line-max-likelihood.png")

# Set up the sampler.
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))

# Clear and run the production chain.
print("Running MCMC...")
sampler.run_mcmc(pos, 500, rstate0=np.random.get_state())
print("Done.")

pl.clf()
fig, axes = pl.subplots(3, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(m_true, color="#888888", lw=2)
axes[0].set_ylabel("$m$")

axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(b_true, color="#888888", lw=2)
axes[1].set_ylabel("$b$")

axes[2].plot(np.exp(sampler.chain[:, :, 2]).T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(f_true, color="#888888", lw=2)
axes[2].set_ylabel("$f$")
axes[2].set_xlabel("step number")

fig.tight_layout(h_pad=0.0)
fig.savefig("line-time.png")

# Make the triangle plot.
burnin = 50
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

fig = triangle.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])
fig.savefig("line-triangle.png")

# Plot some samples onto the data.
pl.figure()
for m, b, lnf in samples[np.random.randint(len(samples), size=100)]:
    pl.plot(xl, m*xl+b, color="k", alpha=0.1)
pl.plot(xl, m_true*xl+b_true, color="r", lw=2, alpha=0.8)
pl.errorbar(x, y, yerr=yerr, fmt=".k")
pl.ylim(-9, 9)
pl.xlabel("$x$")
pl.ylabel("$y$")
pl.tight_layout()
pl.savefig("line-mcmc.png")

# Compute the quantiles.
samples[:, 2] = np.exp(samples[:, 2])
m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""MCMC result:
    m = {0[0]} +{0[1]} -{0[2]} (truth: {1})
    b = {2[0]} +{2[1]} -{2[2]} (truth: {3})
    f = {4[0]} +{4[1]} -{4[2]} (truth: {5})
""".format(m_mcmc, m_true, b_mcmc, b_true, f_mcmc, f_true))






# fit_surface_relation_FID808()
# check_Qmodel_surface_validity_FID808()
check_Qmodel_surface_validity(bolo_name)