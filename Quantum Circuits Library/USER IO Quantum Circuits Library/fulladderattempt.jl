using ITensors

let
    
    print("\nInput 1st Bit: ")
    num1 = readline()
    num1 = parse(Int64, num1) 
    print("Input 2nd Bit: ")
    num2 = readline()
    num2 = parse(Int64, num2) 
    print("Input Carry-In Bit: ")
    num3 = readline()
    num3 = parse(Int64, num3) 
    N = 5
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
        
        psi = MPS(s,"Up")
        v0 = expect(psi, "Proj0"; sites = 1)
        v1 = expect(psi, "Proj1"; sites = 1)
        println("S1: $v0 |0> + $v1 |1>")
  
        v0 = expect(psi, "Proj0"; sites = 2)
        v1 = expect(psi, "Proj1"; sites = 2)
        println("S2: $v0 |0> + $v1 |1>")

        v0 = expect(psi, "Proj0"; sites = 3)
        v1 = expect(psi, "Proj1"; sites = 3)
        println("S2: $v0 |0> + $v1 |1>")

        if num1 == 1
            hj = op("X",s1)
            push!(gates, hj)
        end
        if num2 == 1
            hj = op("X",s2)
            push!(gates, hj)
        end
        if num3 == 1
            hj = op("X",s3)
            push!(gates, hj)
        end
        hj = op("CX",[s1,s4])
        push!(gates, hj)
        hj = op("CX",[s2,s4])
        push!(gates, hj)
        hj = op("CX",[s3,s4])
        push!(gates, hj)
        hj = op("Toffoli",[s1,s2,s5])
        push!(gates, hj)
        hj = op("Toffoli",[s1,s3,s5])
        push!(gates, hj)
        hj = op("Toffoli",[s2,s3,s5])
        push!(gates, hj)
  
        psi = apply(gates, psi; cutoff)
        v0 = expect(psi, "Proj0"; sites = 4)
        v1 = expect(psi, "Proj1"; sites = 4)
        println("Value: $v0 |0> + $v1 |1>")
  
        v0 = expect(psi, "Proj0"; sites = 5)
        v1 = expect(psi, "Proj1"; sites = 5)
        println("Carry: $v0 |0> + $v1 |1>")
        
        v = expect(psi, "Proj1", sites = 4)
        println("Value: $v")
        v = expect(psi, "Proj1", sites = 5)
        println("Carry: $v")
        
    return
  end
  
  