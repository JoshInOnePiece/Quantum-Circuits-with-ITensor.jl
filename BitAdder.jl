using ITensors
using Printf

let
    
    
    print("Input 1st N-Bit Binary Number: ")
    bit1String = readline()
    bit1Array = Array{Int64}(undef, length(bit1String))
    for i = length(bit1String):-1:1
        bit1Array[i] = parse(Int64, bit1String[i])
    end
    
    print("Input 2nd N-Bit Binary Number: ")
    bit2String = readline()
    bit2Array = Array{Int64}(undef, length(bit2String))
    for i = length(bit2String):-1:1
        bit2Array[i] = parse(Int64, bit2String[i])
    end
    
    if length(bit1String) > length(bit2String)
        for i in 1:(length(bit1String)-length(bit2String))
            pushfirst!(bit2Array,0)
        end
    elseif length(bit2String) > length(bit1String)
        for i in 1:(length(bit2String)-length(bit1String))
            pushfirst!(bit1Array,0)
        end
    end

    lengthOfEachBit = length(bit1Array)

    N = lengthOfEachBit*4
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0
    Rz = Array{Float64}(undef, N)
    s = siteinds("Qubit", N)
    gates = ITensor[]
        
        psi = MPS(s,"Up")
        
        #Initialzing Sites
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
        
        psi = apply(gates, psi; cutoff)
        
        result = Array{Int64}(undef, lengthOfEachBit+1)
        print("Answer: ");
        measurement = 4*lengthOfEachBit;
        result[1] = expect(psi, "Proj1", sites = measurement)
        measurement = measurement-1
        result[2] = expect(psi, "Proj1", sites = measurement)
        for i in 3:(lengthOfEachBit+1)
            measurementSum = measurement - 2(i-2)
            result[i] = expect(psi, "Proj1", sites = measurementSum)
        end
        for i in 1:(lengthOfEachBit+1)
            
            @printf("%i", result[i]);
        end
        
    return
  end