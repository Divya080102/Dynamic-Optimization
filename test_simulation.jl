using MyReactor
# define a profile
profile = MyReactor.PiecewiseLinearProfile([0.5 1.0])  # half the time in first interval
# simulate the process
results = MyReactor.simulation(profile)
PyPlot.plot(results.t, [results.u[i][3] for i=1:length(results.u)])
PyPlot.xlabel("t")
PyPlot.ylabel("Product Concentration")
