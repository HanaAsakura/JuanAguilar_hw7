n_steps=20000
n_burn=15000
all: mcmc_lotka_volterra_graphs.pdf

mcmc_lotka_volterra_graphs.pdf:mcmc_lotka_volterra.py
	python mcmc_lotka_volterra.py lotka_volterra_$(n_steps)_$(n_burn).dat &

mcmc_lotka_volterra_$(n_steps)_$(n_burn).dat:mcmc_lotkavolterra.x
	./mcmc_lotkavolterra.x $(n_steps) $(n_burn)

mcmc_lotkavolterra.x:mcmc_lotkavolterra.c
	cc mcmc_lotkavolterra.c -lgsl -lgslcblas -lm -o mcmc_lotkavolterra.x

