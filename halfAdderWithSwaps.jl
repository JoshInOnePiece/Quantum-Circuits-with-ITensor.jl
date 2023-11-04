using ITensors
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
    print("Input 1st Number: ")
    num1 = readline()
    num1 = parse(Int64, num1) 
    print("Input 2nd Number: ")
    num2 = readline()
    num2 = parse(Int64, num2) 
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
        
  
        psi = MPS(s,"Up")
        v0 = expect(psi, "Proj0"; sites = 1)
        v1 = expect(psi, "Proj1"; sites = 1)
        println("S1: $v0 |0> + $v1 |1>")
  
        v0 = expect(psi, "Proj0"; sites = 2)
        v1 = expect(psi, "Proj1"; sites = 2)
        println("S2: $v0 |0> + $v1 |1>")

        if num1 == 1
            hj = op("X",s1)
            push!(gates, hj)
        end
        if num2 == 1
            hj = op("X",s2)
            push!(gates, hj)
        end
        swap(s, gates, 1, 2)
        hj = op("CX",[s[2],s[3]])
        push!(gates, hj)
        swap(s, gates, 1, 2)
        hj = op("CX",[s[2],s[3]])
        push!(gates, hj)
        swap(s, gates, 2, 3)
        swap(s, gates, 1, 2)
        hj = op("Toffoli",[s[2],s[3],s[4]])
        push!(gates, hj)
        swap(s, gates, 1, 2)
        swap(s, gates, 2, 3)
  
        psi = apply(gates, psi; cutoff)
        v0 = expect(psi, "Proj0"; sites = 3)
        v1 = expect(psi, "Proj1"; sites = 3)
        println("Value: $v0 |0> + $v1 |1>")
  
        v0 = expect(psi, "Proj0"; sites = 4)
        v1 = expect(psi, "Proj1"; sites = 4)
        println("Carry: $v0 |0> + $v1 |1>")
        
  
        v = expect(psi, "Proj1", sites = 3)
        println("Value: $v")
        v = expect(psi, "Proj1", sites = 4)
        println("Carry: $v")
        
    return
  end
  
  