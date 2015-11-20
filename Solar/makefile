n_steps=20000
n_burn=15000
all: mcmc_solar_graphs.pdf

mcmc_solar_graphs.pdf:mcmc_solar.py
	python mcmc_solar.py mcmc_solar_$(n_steps)_$(n_burn).dat &

mcmc_solar_$(n_steps)_$(n_burn).dat:mcmc_solar.x
	./mcmc_solar.x $(n_steps) $(n_burn)

mcmc_solar.x:mcmc_solar.c
	cc mcmc_solar.c -lgsl -lgslcblas -lm -o mcmc_solar.x

