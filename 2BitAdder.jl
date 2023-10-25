using ITensors
using Printf

let
    
    print("\nInput 1st 2-Bit Binary Number: ")
    threeBit1 = readline()
    threeBit1Byte1 = parse(Int64,threeBit1[2]);
    threeBit1Byte2 = parse(Int64,threeBit1[1]);
    
    
    print("Input 2nd 2-Bit Binary Number: ")
    threeBit2 = readline()
    threeBit2Byte1 = parse(Int64, threeBit2[2]);
    threeBit2Byte2 = parse(Int64, threeBit2[1]);
    
    N = 8
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0
    Rz = Array{Float64}(undef, N)
    s = siteinds("Qubit", N)
    gates = ITensor[]
        s1 = s[1]
        s2 = s[2]
        s3 = s[3]
        s4 = s[4]
        s5 = s[5]
        s6 = s[6]
        s7 = s[7]
        s8 = s[8]
        
  
        psi = MPS(s,"Up")
        if threeBit1Byte1 == 1
            hj = op("X",s1)
            push!(gates, hj)
        end
        if threeBit1Byte2 == 1
            hj = op("X",s2)
            push!(gates, hj)
        end
        if threeBit2Byte1 == 1
            hj = op("X",s3)
            push!(gates, hj)
        end
        if threeBit2Byte2 == 1
            hj = op("X",s4)
            push!(gates, hj)
        end
        hj = op("CX",[s1,s5])
        push!(gates, hj)
        hj = op("CX",[s3,s5])
        push!(gates, hj)
        hj = op("Toffoli",[s1,s3,s6])
        push!(gates, hj)
        hj = op("CX",[s2,s7])
        push!(gates, hj)
        hj = op("CX",[s4,s7])
        push!(gates, hj)
        hj = op("CX",[s6,s7])
        push!(gates, hj)
        hj = op("Toffoli",[s2,s4,s8])
        push!(gates, hj)
        hj = op("Toffoli",[s2,s6,s8])
        push!(gates, hj)
        hj = op("Toffoli",[s4,s6,s8])
        push!(gates, hj)
  
        psi = apply(gates, psi; cutoff)
        
        LSB = expect(psi, "Proj1", sites = 5)
        MSB = expect(psi, "Proj1", sites = 7)
        Carry = expect(psi, "Proj1", sites = 8)
        println("Carry: $Carry")
        println("MSB: $MSB")
        println("LSB: $LSB")
        @printf("Answer: %i%i%i", Carry, MSB, LSB)
    return
  end