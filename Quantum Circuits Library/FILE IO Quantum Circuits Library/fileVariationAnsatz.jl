using ITensors
using Plots
using Statistics
function fileVariationalCircuit(inputFile, outputFile)
    include(raw"Functions Library\File IO Functions\fileEntanglementFunctions.jl")
    let
        #Reads inputs values
        numOfSteps = readline(inputFile)
        timeSteps = readline(inputFile)
        totalTime = readline(inputFile)

        #Parses input to integers and floats
        numOfStepsInt = parse(Int64,numOfSteps)
        timeStepsFloat = parse(Float64,timeSteps)
        totalTimeFloat = parse(Float64,totalTime)

        N = numOfStepsInt
        cutoff = 1E-8
        tau = timeStepsFloat
        ttotal = totalTimeFloat
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
            measurementArray = Array{Float64}(undef, Int.(ttotal/tau)+1)
            timeStepsArray = Array{Float64}(undef, Int.(ttotal/tau)+1)
            counter = 1
            # Compute and print <Sz> at each time step
            # then apply the gates to go to the next time
            for t in 0.0:tau:ttotal
                Rz = Array{Float64}(undef, N)
                for c in 1:1:N
                    Rz[c] = expect(psi, "Sz"; sites=c)
                    average = mean(Rz)
                    write(outputFile, "Site $c at $t: $average \n")
                    if(c == N)
                        timeStepsArray[counter] = t
                        measurementArray[counter] = average
                        counter = counter + 1
                    end
                    #println("Site $c at $t: $average")
                end
                t≈ttotal && break
                psi = apply(gates, psi; cutoff)
                normalize!(psi)
            end
            filePromptForMeasuringEntanglement(inputFile, outputFile, N, psi)
            scatter(timeStepsArray, measurementArray, xlabel="Time", ylabel="S_z", show = true, label = "Sz", title = "Sz values of Variational Ansantz Circuit over Time")
            #If you want a linear line
            #plot(timeStepsArray, measurementArray, xlabel="Time", ylabel="S_z", show = true) if yo
            savefig(raw"IO FILES\Scatterplots\VariationalAnsatzScatterPlot.png")
        return   
    end
end