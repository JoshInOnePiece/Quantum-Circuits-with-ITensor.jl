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