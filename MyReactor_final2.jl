module MyReactor

using DifferentialEquations     # for solving the system of ODEs
using Fresa                     # for optimization 
using Printf                    # for formatted output
using PyPlot                    # for plotting profiles

const tfinal = 120.0              # time in hours
const Vmin = 60.0              # in Litres
const Vmax = 100.0              # in Litres

abstract type VolumeProfile end

mutable struct PiecewiseLinearProfile <: VolumeProfile	#mutable structure to store feedrate vs time profile
    ft :: Array{Float64}        # time deltas
    F :: Array{Float64}        # temperature deltas
    V0 :: Float64               # initial temperature
    function PiecewiseLinearProfile(ft, F, V0)
        if length(ft) != length(F)		#check for equal number of feedrate and time intervals
            error("Number of time points and temperature changes must match.")
        end
        if V0 < Vmin || V0 > Vmax
            error("V0 = $(V0) not in [$Vmin,$Vmax]")
        end
        new(ft,F,V0)
    end
    function PiecewiseLinearProfile(t :: PiecewiseLinearProfile)
        new(t.ft, t.F, t.V0)
    end
end
# allow up to 10 degrees change per time interval
ΔV_max = 40.0

function Feedrate(p:: PiecewiseLinearProfile, tf :: Float64)		#function to return feedrate at any given time
     F = 0.0
     flag = 0
     for i=1:(length(p.ft)-1)
	if(tf>=p.ft[i] && tf<p.ft[i+1])
		F = p.F[i]
		flag = 1
	end
     end
     if(flag == 0)
	F = p.F[length(p.ft)]
     end
     return F
end
    

function Vl(profile :: PiecewiseLinearProfile, t :: Float64) :: Float64		#function to return volume at any given time using the feedrate 
    τ = 0 
    V = profile.V0
    for i=1:length(profile.ft)
	if(i == length(profile.ft))
		Δτ = tfinal - profile.ft[i]
	else
        	Δτ = profile.ft[i+1] - profile.ft[i] 
	end
	
        ΔV = profile.F[i] * Δτ
#	println("delta t, V = $Δτ $ΔV")
        if t < τ + Δτ
            # partial step to reach desired time
            V += ΔV * (t-τ)/Δτ
            # correct T if out of bounds
            if V < Vmin
                V = Vmin
            elseif V>Vmax
                V = Vmax
            end
            # leave loop
            break
        else
            # take full step
            V += ΔV
            # correct V if out of bounds
            if V < Vmin
                V = Vmin
            elseif V>Vmax
                V = Vmax
            end
        end
        τ += Δτ
    end
   # println("$V $t")
    return V
end

#function to formulate the mathematical equations of the paper
function reactor(delements, elements :: Array{Float64}, 	#the array elements contains values of X(microbial growth concentration), P(concentration of Product),S(concentration of substrate), CO(dissolved oxygen content) respectively
		 p :: PiecewiseLinearProfile,			
		 tf:: Float64
                 ) 
     F = Feedrate(p,tf)			
     V = Vl(p,tf)

    		#constant CO concentration in the reactor
    # println("$t")
    X, P, S, CO = elements
    
    delements[1] = ((5.53537543*(10^-12)*(S/(S + 0.1828*X))*(CO/(0.0352*X + CO))*(1-(X/0.87))) - (4.60796599*(10^-62)*(1-(CO/(0.0368+CO)))) - F/V)*X
    delements[2] =  0.05*S*X/(0.0002+S+(10*S^2)) - 4*10^(-4)*P - F*P/V	
    delements[3] = -0.0624*X - 4*delements[1] - ((1/0.68)*delements[2]) - F*X/V 
    delements[4] =  ((1.2720140414/(V^0.4))*(0.037 - CO)) - (0.004*X) - ((1/43.5)*delements[1]) - ((1/253.3)*delements[2]) - (F*CO/V)	#created the ODE equations
end

function Vplot(profile :: PiecewiseLinearProfile)	#function to plot volume profile
    n = 120
   # δx = 1.0/1000
    x = [1.0*i for i in 1:n+1]
    y = [1.0*Vl(profile,x[i]) for i in 1:n+1]

    PyPlot.plot(x,y, linewidth=2.0, linestyle="--", )
    PyPlot.ylabel("Volume (K)")
    PyPlot.xlabel("Time")
    PyPlot.title("Volume profile")
end

function Fplot(profile :: PiecewiseLinearProfile)	#function to plot feedrate profile
    n = 120
   # δx = 1.0/1000
    x = [2.0*i for i in 1:n+1]
    y = [1.0*Feedrate(profile,x[i]) for i in 1:n+1]

    PyPlot.plot(x,y)
    PyPlot.ylabel("FeedRate (K)")
    PyPlot.xlabel("Time")
    PyPlot.title("Feedrate profile")
   # PyPlot.xlim(0,200)
end


function simulation(profile :: PiecewiseLinearProfile)		#do the simulation
    # for testing some performance aspects, we have two ways of
    # simulating the reactor, i.e. solving the differential equations:
    # an efficient solver (RK 2-3) and a computationally inefficient
    # one, Euler's method
    efficient = true
    if efficient
        # efficient solver
        tspan = (0.0,tfinal)
	 	
        prob = ODEProblem(reactor, [0.05, 0.0, 10.0, 0.01], tspan, profile)
        results = DifferentialEquations.solve(prob,saveat = 19.5)
        results
    else
        δt = tfinal/1e6
        t = 0.0
        var :: Array{Float64} = [0.05, 0.0, 10.0, 0.01]
        while t < tfinal
            if t+δt > tfinal
                δt = tfinal - t
            end
	 #   V = Vl(profile,t+δt)
	 #   if(V>Vmax)
#		break
#	    end
            δvar = δt * reactor(var,profile,t+δt)
            var = var + δvar
	    t = t + δt
        end
        var                       # end values
    end
end

end

