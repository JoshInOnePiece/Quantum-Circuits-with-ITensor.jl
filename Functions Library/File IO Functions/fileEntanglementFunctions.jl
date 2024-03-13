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
function printEverySiteEntanglement(N,psi,outputFile)
    for i in 2:N
        entanglement = calcEntanglement(i,psi)
        write(outputFile, "\nSite $i: $entanglement")
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
function filePromptForMeasuringEntanglement(inputFile, outputFile, N, psi)
        write(outputFile,"\n\n")
		numSiteString = readline(inputFile)
		numSite = parse(Int64, numSiteString)
		entanglementOfSetOfSites = calcEntanglement(numSite, psi)
		write(outputFile,"Entanglement between sites 1 - $numSite: $entanglementOfSetOfSites\n")
		write(outputFile,"Total Entanglement: ")
		totalEntanglement = summationOfSiteEntanglement(N, psi)
		write(outputFile,"$totalEntanglement\n")
		write(outputFile,"Average Entanglement: ")
		averageEntanglement = totalEntanglement/N
		write(outputFile,"$averageEntanglement\n")
end