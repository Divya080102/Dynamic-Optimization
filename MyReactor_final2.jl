module MyReactor

using DifferentialEquations     # for solving the system of ODEs
using Fresa                     # for optimization 
using Printf                    # for formatted output
using PyPlot                    # for plotting profiles

#constant conditions of the reactor - 
const tfinal = 390.0		# hours
const T = 303.0			# Kelvin
const pH = 7.0			#Value of pH for the reactions
const F = 0.5			#feedrate
const Vmax = 100.0		#Maximum volume of contents of the reactor(= volume of the reactor container)
const Vmin = 60.0		#Initial volume of the reactor

mutable struct PiecewiseLinearProfile	#Structure to create profile for time spent in the batch reactor
    ft :: Array{Float64}        # time deltas
    function PiecewiseLinearProfile(ft)	
	new(ft)
    end
    function PiecewiseLinearProfile(t :: PiecewiseLinearProfile)
        new(t.ft)
    end
end

function volume(profile :: PiecewiseLinearProfile, t :: Float64) :: Float64	#Finding volume at a particular time point
										#Rate of volume change = Feedrate(F) = 0.5
    τ = 0
    V = Vmin	
    for i=1:length(profile.ft)
        Δτ = profile.ft[i] * (tfinal - τ)
	ΔV = F*Δτ
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
            # correct T if out of bounds
            if V < Vmin
                V = Vmin
            elseif V>Vmax
                V = Vmax
            end
        end
        τ += Δτ
    end
    return V
end

#function to formulate the mathematical equations of the paper
function reactor(delements, elements :: Array{Float64}, 	#the array elements contains values of X(microbial growth concentration), S(concentration of substrate), P(concentration of Product) respectively
		 p :: PiecewiseLinearProfile,			
		 tf:: Float64
                 ) 
   # t = time(p,tf)
    V = volume(p,tf)		#finding the volume of reactor contents at a given time tf
 #   println("Volume = $V")
    CO = 0.037		#constant CO concentration in the reactor
    # println("$t")
    X, S, P = elements
    delements[1] = ((5.53537543*(10^-12)*(S/(S + 0.1828*X))*(CO/(0.0352*X + CO))*(1-(X/0.87))) - (4.60796599*(10^-62)*(1-(CO/(0.0368+CO)))) - F/V)*X
    delements[2] = -0.0624*X - 4*(((5.53537543*(10^-12)*(S/(S + 0.1828*X))*(CO/(0.0352*X + CO))*(1-(X/0.87))) - (4.60796599*(10^-62)*(1-(CO/(0.0368+CO)))) - F/V)*X) - ((1/0.68)* 0.05*S*X/(0.0002+S+(10*S^2)) - 4*10^(-4)*P - F*P/V) - F*X/V
    delements[3] =  0.05*S*X/(0.0002+S+(10*S^2)) - 4*10^(-4)*P - F*P/V			#created the ODE equations 
  #  delements[4] =  F	
 #   print("X = $X, S = $S, P = $P")
   # return delements
end



function simulation(profile :: PiecewiseLinearProfile)
    # for testing some performance aspects, we have two ways of
    # simulating the reactor, i.e. solving the differential equations:
    # an efficient solver (RK 2-3) and a computationally inefficient
    # one, Euler's method
    efficient = true
    if efficient
        # efficient solver
        tspan = (0.0,tfinal)
        prob = ODEProblem(reactor, [0.05, 40.0, 0.0], tspan, profile)
        results = DifferentialEquations.solve(prob,saveat = 19.5)
        results
    else
        δt = tfinal/1e4
        t = 0.0
        var :: Array{Float64} = [0.05, 40.0, 0.0]
        while t < tfinal
            if t+δt > tfinal
                δt = tfinal - t
            end
            δvar = δt * reactor(var,profile,t+δt)
            var = var + δvar
        end
        var                       # end values
    end
end

end
