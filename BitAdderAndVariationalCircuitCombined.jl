using ITensors
using Statistics
using Printf

function adderCircuit(inputFile, outputFile)
    let
        #reads in input
        bit1String = " "
        bit2String = " "
        bit1String = readline(inputFile)
        bit2String = readline(inputFile)
        
        #Parses input into and integer and stores it into an array
        bit1Array = Array{Int64}(undef, length(bit1String))
        for i = length(bit1String):-1:1
            bit1Array[i] = parse(Int64, bit1String[i])
        end
        
        #Same as above but for the second binary bits
        bit2Array = Array{Int64}(undef, length(bit2String))
        for i = length(bit2String):-1:1
            bit2Array[i] = parse(Int64, bit2String[i])
        end
        
        #Checks the sizes of both arrays, and makes the smaller one the same size as the larger by padding it with zeroes
        if length(bit1String) > length(bit2String)
            for i in 1:(length(bit1String)-length(bit2String))
                pushfirst!(bit2Array,0)
            end
        elseif length(bit2String) > length(bit1String)
            for i in 1:(length(bit2String)-length(bit1String))
                pushfirst!(bit1Array,0)
            end
        end

        #Length of each bit
        lengthOfEachBit = length(bit1Array)

        #Initialzing site indices
        N = lengthOfEachBit*4
        cutoff = 1E-8
        tau = 0.1
        ttotal = 5.0
        Rz = Array{Float64}(undef, N)
        s = siteinds("Qubit", N)
        gates = ITensor[]
            
            psi = MPS(s,"Up")
            
            #Initialzing Sites based on array
            for i = (lengthOfEachBit):-1:1
                if bit1Array[i] == 1
                    index = lengthOfEachBit-(i-1)
                    hj = op("X",s[index])
                    push!(gates, hj)
                end
                
            end
            
            for i = (lengthOfEachBit):-1:1
                if bit2Array[i] == 1
                    index = lengthOfEachBit+lengthOfEachBit-(i-1)
                    hj = op("X",s[index])
                    push!(gates, hj)
                end
            end
            
            #Half Adder Code
            hj = op("CX",[s[1],s[lengthOfEachBit*2+1]])
            push!(gates, hj)
            hj = op("CX",[s[lengthOfEachBit+1],s[lengthOfEachBit*2+1]])
            push!(gates, hj)
            hj = op("Toffoli",[s[1],s[lengthOfEachBit+1],s[2*lengthOfEachBit+2]])
            push!(gates, hj)

            #Full Adder Code
            for i in 1:(lengthOfEachBit-1)

                hj = op("CX",[s[1+i],s[2*lengthOfEachBit+2i+1]])
                push!(gates, hj) 
                hj = op("CX",[s[lengthOfEachBit+i+1],s[2*lengthOfEachBit+2*i+1]])
                push!(gates, hj) 
                hj = op("CX",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1]])
                push!(gates, hj)
                hj = op("Toffoli",[s[1+i],s[1+i+lengthOfEachBit],s[2*lengthOfEachBit+2*i+2]])
                push!(gates, hj)
                hj = op("Toffoli",[s[1+i],s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+2]])
                push!(gates, hj)
                hj = op("Toffoli",[s[1+i+lengthOfEachBit],s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+2]])
                push!(gates, hj) 
            end
            
            #Opens file where output will be stored 
            psi = apply(gates, psi; cutoff)

            #Result Array will be used to store the output in the file and print it on the terminal
            result = Array{Int64}(undef, lengthOfEachBit+1)
            print("Answer: ");

            #Iterating through the sites and measuring them to get the output
            measurement = 4*lengthOfEachBit;
            result[1] = expect(psi, "Proj1", sites = measurement)
            measurement = measurement-1
            result[2] = expect(psi, "Proj1", sites = measurement)
            for i in 3:(lengthOfEachBit+1)
                measurementSum = measurement - 2(i-2)
                result[i] = expect(psi, "Proj1", sites = measurementSum)
            end

            #Outputs sum to the output and prints it on the terminal
            for i in 1:(lengthOfEachBit+1)
                write(outputFile, string(result[i]))
                @printf("%i", result[i]);
            end
            println("")
    return
    end
end
function variationalCircuit(inputFile, outputFile)
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
            # Compute and print <Sz> at each time step
            # then apply the gates to go to the next time
            for t in 0.0:tau:ttotal
                for c in 1:1:N
                    Rz[c] = expect(psi, "Sz"; sites=c)
                    average = mean(Rz)
                    write(outputFile, "Site $c at $t: $average \n")
                    println("Site $c at $t: $average")
                end
                t≈ttotal && break
                psi = apply(gates, psi; cutoff)
                normalize!(psi)
            end
        return   
    end
end
let
    #Names of input and output files
    inputAdderFileName = raw"C:\Users\manik\OneDrive\Documents\Josh's Quantum World\adderInput.txt"
    outputAdderFileName = raw"C:\Users\manik\OneDrive\Documents\Josh's Quantum World\adderOutput.txt"
    inputVariationalFileName = raw"C:\Users\manik\OneDrive\Documents\Josh's Quantum World\variationalInput.txt"
    outputVariationalFileName = raw"C:\Users\manik\OneDrive\Documents\Josh's Quantum World\variationalOutput.txt"

    #opens files from which input and output is stored
    inputAdderFile = open(inputAdderFileName,"r")
    outputAdderFile = open(outputAdderFileName, "w")
    inputVariationalFile = open(inputVariationalFileName,"r")
    outputVariationalFile = open(outputVariationalFileName, "w")

    #Calling the Binary Adder Circuit
    adderCircuit(inputAdderFile, outputAdderFile)
    variationalCircuit(inputVariationalFile, outputVariationalFile)
    
    #Closes files
    close(inputAdderFile)
    close(inputVariationalFile)
    close(outputAdderFile)
    close(outputVariationalFile)
end