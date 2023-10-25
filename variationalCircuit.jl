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
			println("$t $average")
		end
		t≈ttotal && break
		psi = apply(gates, psi; cutoff)
		normalize!(psi)
end

return
end