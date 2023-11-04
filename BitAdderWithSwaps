#=Similar to the other binary adder except that the circuit has the restriction
where only adjacent site indices can be connected by a gate.
In order to connect gates of sites that are not adjacent is to swap the value with another site
=#
using ITensors
using Printf
function swap(s, gates, start, ending)
    let 
        psi = MPS(s,"Up")
        for i in (start):(ending-1)
            hj = op("SWAP",[s[i],s[i+1]])
            push!(gates, hj)
        end 
    end
end

function unswap(s, gates, start, ending)
    let 
        psi = MPS(s,"Up")
        for i = (ending):-1:(start+1)
            hj = op("SWAP",[s[i],s[i-1]])
            push!(gates, hj)
        end
        
    end
end
let
    #takes in input for first binary number and puts it into an array
    print("Input 1st N-Bit Binary Number: ")
    bit1String = readline()
    bit1Array = Array{Int64}(undef, length(bit1String))
    for i = length(bit1String):-1:1
        bit1Array[i] = parse(Int64, bit1String[i])
    end

    #takes input for second binary number and puts it into an array
    print("Input 2nd N-Bit Binary Number: ")
    bit2String = readline()
    bit2Array = Array{Int64}(undef, length(bit2String))
    for i = length(bit2String):-1:1
        bit2Array[i] = parse(Int64, bit2String[i])
    end

    #Finds out which string is longer and sets the array of the smaller number to the same size as the larger number
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
        swap(s, gates, 1, lengthOfEachBit*2)
        hj = op("CX",[s[lengthOfEachBit*2],s[lengthOfEachBit*2+1]])
        push!(gates, hj)
        unswap(s, gates, 1, lengthOfEachBit*2)

        swap(s, gates, lengthOfEachBit+1, lengthOfEachBit*2)
        hj = op("CX",[s[2*lengthOfEachBit],s[lengthOfEachBit*2+1]])
        push!(gates, hj)
        unswap(s, gates, lengthOfEachBit+1, lengthOfEachBit*2)

        swap(s, gates, lengthOfEachBit+1, 2*lengthOfEachBit+1)
        swap(s, gates, 1, 2*lengthOfEachBit)
        hj = op("Toffoli",[s[2*lengthOfEachBit],s[2*lengthOfEachBit+1],s[2*lengthOfEachBit+2]])
        push!(gates, hj)
        unswap(s, gates, 1, 2*lengthOfEachBit)
        unswap(s, gates, lengthOfEachBit+1, 2*lengthOfEachBit+1)

        #Full adder with swaps
        for i in 1:(lengthOfEachBit-1)
            swap(s, gates, 1+i, 2*lengthOfEachBit+2i)
            hj = op("CX",[s[2*lengthOfEachBit+2i],s[2*lengthOfEachBit+2i+1]])
            push!(gates, hj) 
            unswap(s, gates, 1+i, 2*lengthOfEachBit+2i)

            swap(s, gates, lengthOfEachBit+i+1, 2*lengthOfEachBit+2*i)
            hj = op("CX",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1]])
            push!(gates, hj)
            unswap(s, gates, lengthOfEachBit+i+1, 2*lengthOfEachBit+2*i)
            
            hj = op("CX",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1]])
            push!(gates, hj)

            swap(s, gates, 1+i+lengthOfEachBit, 2*lengthOfEachBit+2*i+1)
            swap(s, gates, 1+i, 2*lengthOfEachBit+2*i)
            hj = op("Toffoli",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1],s[2*lengthOfEachBit+2*i+2]])
            push!(gates, hj)
            unswap(s, gates, 1+i, 2*lengthOfEachBit+2*i)
            unswap(s, gates, 1+i+lengthOfEachBit, 2*lengthOfEachBit+2*i+1)

            swap(s, gates, 2*lengthOfEachBit+2*i, 2*lengthOfEachBit+2*i+1)
            swap(s, gates, 1+i, 2*lengthOfEachBit+2*i)
            hj = op("Toffoli",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1],s[2*lengthOfEachBit+2*i+2]])
            push!(gates, hj)
            unswap(s, gates, 1+i, 2*lengthOfEachBit+2*i)
            unswap(s, gates, 2*lengthOfEachBit+2*i, 2*lengthOfEachBit+2*i+1)

            swap(s, gates, 2*lengthOfEachBit+2*i, 2*lengthOfEachBit+2*i+1)
            swap(s, gates, 1+i+lengthOfEachBit, 2*lengthOfEachBit+2*i)
            hj = op("Toffoli",[s[2*lengthOfEachBit+2*i],s[2*lengthOfEachBit+2*i+1],s[2*lengthOfEachBit+2*i+2]])
            push!(gates, hj) 
            unswap(s, gates, 1+i+lengthOfEachBit, 2*lengthOfEachBit+2*i)
            unswap(s, gates, 2*lengthOfEachBit+2*i, 2*lengthOfEachBit+2*i+1)
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
