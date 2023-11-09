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
        #print("\nSite $i: $total()")
    end
    #print("\nTotal Entanglement: $total")
    return total
end