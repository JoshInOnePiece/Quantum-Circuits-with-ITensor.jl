using ITensors
using Plots
using Statistics
function fileTEBDCircuit(inputFile, outputFile)
  include(raw"Functions Library\File IO Functions\fileEntanglementFunctions.jl")
  let
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
    s = siteinds("S=1/2", N; conserve_qns=true)
  
    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Initialize psi to be a product state (alternating up and down)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  
    c = div(N, 2) # center site
    measurementArray = Array{Float64}(undef, Int.(ttotal/tau)+1)
    timeStepsArray = Array{Float64}(undef, Int.(ttotal/tau)+1)
    counter = 1
    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    write(outputFile, "Sz measurement at each time step: \n")
    for t in 0.0:tau:ttotal
      Sz = expect(psi, "Sz"; sites=c)
      write(outputFile,"$t: $Sz\n")
      measurementArray[counter] = Sz
      timeStepsArray[counter] = t
      counter = counter + 1
      t≈ttotal && break
      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    #Measures the Entanglement across a certain number of sites
    filePromptForMeasuringEntanglement(inputFile,outputFile, N, psi)
    scatter(timeStepsArray, measurementArray, xlabel="Time", ylabel="S_z", show = true, label = "Sz", title = "Sz values of TEBD Circuit over Time")
    #If you want a line
    #plot(timeStepsArray, measurementArray, xlabel="Time", ylabel="S_z", show = true)
    #Saves plot into this file
    savefig(raw"IO FILES\Scatterplots\TEBDScatterPlot.png")
    return
  end
end