using ITensors

function calcEntanglement(siteIndex, psi)
    let
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
end
function printEverySiteEntanglement(N,psi)
    let 
        for i in 2:N
            entanglement = calcEntanglement(i,psi)
            print("\nSite $i: $entanglement")
        end
    end
end
function summationOfSiteEntanglement(N, psi)
    let
        total = 0.0;
        for i in 2:N
            total += calcEntanglement(i, psi)
            #print("\nSite $i: $total()")
        end
        #print("\nTotal Entanglement: $total")
        return total
    end
end
let
  N = 100
  cutoff = 1E-8
  tau = 0.1
  ttotal = 5.0

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

  # Compute and print <Sz> at each time step
  # then apply the gates to go to the next time
  for t in 0.0:tau:ttotal
    Sz = expect(psi, "Sz"; sites=c)
    println("$t $Sz")

    t≈ttotal && break

    psi = apply(gates, psi; cutoff)
    normalize!(psi)
  end

  print("\n\n")
  print("What entanglement: ")
  numSiteString = readline()
  numSite = parse(Int64, numSiteString)
  entanglementOfSetOfSites = calcEntanglement(numSite, psi)
  print("Entanglement between sites 1 - $numSite: $entanglementOfSetOfSites\n")
  print("Total Entanglement: ")
  totalEntanglement = summationOfSiteEntanglement(N, psi)
  print("$totalEntanglement\n")
  print("Average Entanglement: ")
  averageEntanglement = totalEntanglement/N
  print("$averageEntanglement\n")
  return
end