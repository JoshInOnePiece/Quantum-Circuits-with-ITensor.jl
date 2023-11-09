using ITensors
using Statistics
function calcEntanglement(siteIndex, psi)
    entanglement = 0.0
	orthogonalize!(psi, siteIndex)
    U,S,V = svd(psi[siteIndex], (linkind(psi, siteIndex-1), siteind(psi,siteIndex)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        entanglement -= p * log(p)
    end

    return entanglement
end
function printEverySiteEntanglement(N,psi)
    for i in 2:N
        entanglement = calcEntanglement(i,psi)
        print("\nSite $i: $entanglement")
    end
end
function summationOfSiteEntanglement(N, psi)
    total = 0.0;
    for i in 2:N
        total += calcEntanglement(i, psi)
        #print("\nSite $i: $total")
    end
    #print("\nTotal Entanglement: $total")
    return total
end
let
	N = 5
	cutoff = 1E-8
	tau = 0.1
	ttotal = 5.0
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
				println("Site $c at $t: $average")
			end
			t≈ttotal && break
			psi = apply(gates, psi; cutoff)
			normalize!(psi)
		end

		print("\n\n")
        printEverySiteEntanglement(N,psi);
        totalEntanglement = summationOfSiteEntanglement(N, psi);
		averageEntanglement = totalEntanglement/(N-1)
		print("\nTotal Entanglement: $averageEntanglement")
	return
end