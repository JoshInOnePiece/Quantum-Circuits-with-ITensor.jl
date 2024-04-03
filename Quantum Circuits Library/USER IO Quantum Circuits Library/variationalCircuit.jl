using ITensors
using Statistics
using Plots
using StatsPlots
using DataFrames, CSV
using CSV, Tables
using DelimitedFiles
include(raw"Functions Library\Terminal Functions\entanglementFunctions.jl")

let

	#Reads inputs values
    print("Number of Sites: ")
    numOfSteps = readline()
    print("Duration of Time Step: ")
    timeSteps = readline()
    print("Total amount of Time: ")
    totalTime = readline()

	#Parses input to integers and floats
	numOfStepsInt = parse(Int64,numOfSteps)
	timeStepsFloat = parse(Float64,timeSteps)
	totalTimeFloat = parse(Float64,totalTime)

	N = numOfStepsInt #5
	cutoff = 1E-8
	tau = timeStepsFloat #0.1
	ttotal = totalTimeFloat #5.0
	# Make an array of 'site' indices
	s = siteinds("Qubit", N; conserve_qns=false)
	# Make gates (1,2),(2,3),(3,4),...
	gates = ITensor[]
		for j in 1:(N - 1)
			s1 = s[j]
			s2 = s[j + 1]
			hj =
			op("CX",[s1,s2]) +
			op("Rx", s1, θ = rand(Float64)*pi)*op("Rx", s2, θ = rand(Float64)*pi) +
			op("Rz", s1, ϕ = rand(Float64)*pi)*op("Rz", s2, ϕ = rand(Float64)*pi)
			Gj = exp(-im * tau / 2 * hj)
			push!(gates, Gj)
		end

		# Include gates in reverse order too
		# (N,N-1),(N-1,N-2),...
		append!(gates, reverse(gates))

		# Initialize psi to be a product state (alternating up and down)
		psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
		Rz = Array{Float64}(undef, N)
		measurementArray = Array{Any}(undef, Int.(ttotal/tau)+2, 2)
    	measurementArray[1, 1] = "Time"
    	measurementArray[1, 2] = "Measurement"
    	counter = 2
		# Compute and print <Sz> at each time step
		# then apply the gates to go to the next time
		for t in 0.0:tau:ttotal
			for c in 1:1:N
				Rz[c] = expect(psi, "Sz"; sites=c)
				average = mean(Rz)
				println("Site $c at $t: $average")
				if(c == N)
					measurementArray[counter, 1] = Float64.(t)
					measurementArray[counter, 2] = Float64.(average)
					counter = counter + 1
				end
			end
			t≈ttotal && break
			psi = apply(gates, psi; cutoff)
			normalize!(psi)
		end
		promptForMeasuringEntanglement(N, psi)
		outputCSV = raw"IO FILES\CSV FILES\variationalAnsantzOutput.csv"
        writedlm(raw"IO FILES\CSV FILES\variationalAnsantzOutput.csv",  measurementArray, ',')
        data = CSV.read(raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\CSV FILES\variationalAnsantzOutput.csv", DataFrames.DataFrame)
        Plots.scatter(data.Time, data.Measurement, xlabel="Time", ylabel="Measurement", show = true)
	return
end