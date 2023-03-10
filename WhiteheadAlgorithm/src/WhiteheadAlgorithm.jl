module WhiteheadAlgorithm

include("alphabets.jl")
include("freerewrite.jl")
include("automorphisms.jl")

#This function checks if a word w0 is primitve in the free group corresponding to A i.e. if w0 can be part of a basis
#If this is the case we also return an automorphism mapping the first generator to w0
#If pr==true we also print all the intermediate steps
function isPrimitive(A::FreeAlphabet, w0::Vector{Int}, pr=true::Bool)
    if length(freerewrite(A,w0))==0
        println("The trivial element is not primitive")
        return (false, f)
    end
    (v,f)=minimize(A,[w0],pr)
    w=v[1]
    if length(w)==1 #w0 is primitive
        (v,f)=permuteGensandInv(A,v,f,pr) #we may need some extra automorphisms to ensure that the first generator gets mapped to w0
        w=v[1]
        pr&&println(word(A, w0)*" is a primitive word. Here is an automorphism mapping the first generator to it:")
        pr&&print(f)
        return (true,f)
    else
        pr&&println(word(A, w0)*" is not a primitive word, the minimal element in its equivalence class is "*word(A, w))
        return (false, f)
    end
end

#This function checks if the words stored in b0 can be extended to a basis of the free group corresponding to A
#If this is the case we also return an automorphism mapping the first generators to the words in b0
#If pr==true we also print all the intermediate steps
function isPartOfBasis(A::FreeAlphabet, b0::Vector{Vector{Int}}, pr=true::Bool) 
    (b,f)=minimize(A,b0,pr)
    if summedweight(b)==length(b)
        (b,f)=permuteGensandInv(A,b,f,pr) #we may need some extra automorphisms to ensure that the first generators get mapped to b0
        w=v[1]
        lb = length(b)
        pr&&print((w -> word(A,w)).(b0))
        pr&&println(" can be part of a basis. Here's an automorphism mapping the first $lb generators to them:")
        pr&&print(f)
        return (true, f)
    else
        pr&&print((w -> word(A,w)).(b0))
        pr&&println(" can not be part of a basis.")
        return (false,f)
    end
end

#This function takes a vector of words and minimizes their summed lengths using first Nielsen and then Whitehead automorphisms
#If pr==true we also print all the intermediate steps
function minimize(A::FreeAlphabet, b0::Vector{Vector{Int}}, pr=true::Bool)
    b = (w -> freerewrite(A,w)).(b0)
    pr&&println((w -> word(A,w)).(b0))
    pr&&println("")
    n = length(A.gen)
    f = Endomorphism(A, [[A.gen[k]] for k in 1:n]) #f is initialized as Identity
    autfound = true
    for w in b
        if length(w)==0
            println("The trivial element can not be part of a basis")
            return (b, f)
        end
        if sum([w==w2 for w2 in b]) > 1
            println("This cannot be part of a basis because "*word(A,w)*" is contained more than once.")
            return (b, f)
        end
    end
    @assert length(b0) <= n "Too many elements to be a basis of the free group in $n generators"
    #repeat till no length shortening automorphism is found or every word in the vector has length 1
    while autfound && summedweight(b)!=length(b) 
        autfound = false
        #iterate over all Nielsen automs except for x |-> x^-1 because they don't shorten any word
        for (i,j,lr,pm) in Iterators.product(1:n,1:n,['l','r'],['+','-']) 
            if i != j
                NAut = NielsenAut(A, i, j, lr, pm)
                b2 = (w -> evaluate(NAut, w)).(b)
                if summedweight(b2) < summedweight(b) #the automorphism shortened b
                    b = b2
                    pr&&println(NAut)
                    pr&&println((w -> word(A,w)).(b))
                    pr&&println("")
                    pm2 = (pm == '+') ? '-' : '+'
                    #f2 is the inverse of the Nielsen autom used above to shorten b                    
                    f2 = NielsenAut(A, i, j, lr, pm2) 
                    f = compose(f,f2)
                    autfound = true
                    break
                end
            end
        end
        if autfound
            continue
        end
        pr&&println("Now we need to look at all Whitehead automorphisms...")
        for i in 1:n
            if autfound
                break
            end
            #iterate over all Whitehead automorphisms except for the permutations because they don't change the length
            for options in Iterators.product([collect(1:4) for k in 1:(n-1)]...) 
                options = collect(options)
                WAut = WhiteheadAut(A, i, options)
                b2 = (w -> evaluate(WAut, w)).(b)
                invopt = inverseoptions(options)
                WAut2 = WhiteheadAut(A, i, invopt) #WAut2 is the inverse of WAut             
                if summedweight(b2) < summedweight(b)
                    b = b2
                    pr&&println(WAut)
                    pr&&println((w -> word(A,w)).(b))
                    f = compose(f,WAut2)
                    autfound = true
                    break
                end
                #we also need to check whether the inverse WAut2 shortens b
                b3 = (w -> evaluate(WAut2, w)).(b)
                if summedweight(b3) < summedweight(b)
                    b = b3
                    pr&&println(WAut2)
                    pr&&println((w -> word(A,w)).(b))
                    f = compose(f,WAut)
                    autfound = true
                    break
                end
            end
        end
    end
    return(b,f)
end

#We need this function if minimize() managed to minimize the word vector s.t. now every word has length 1
#but some generators are still inverted and their ordering is incorrect
function permuteGensandInv(A::FreeAlphabet, b::Vector{Vector{Int}}, f::Endomorphism, pr=true::Bool)
    pr&&println("")
    #First we map our word vector to another one where each word is a generator 
    #(so no inverses of generators anymore)
    for i in 1:length(b)
        (k, isgen) = A.letterToGen[b[i][1]]
        if !isgen
            Aut = InvertGenAut(A, k)
            b = (w -> evaluate(Aut, w)).(b)
            f = compose(f, Aut)
            pr&&println(Aut)
            pr&&println((w -> word(A,w)).(b))
            pr&&println("")
        end
    end
    #Now we need to get the generators in the right order
    #To achieve that we use transpositions
    for j in 1:length(b)
        k = A.letterToGen[b[j][1]][1]
        if k != j
            T = TranspositionAut(A, k, j)
            b = (w -> evaluate(T, w)).(b)
            f = compose(f, T)
            pr&&println(T)
            pr&&println((w -> word(A,w)).(b))
            pr&&println("")
        end
    end
    return (b,f)
end

end # module WhiteheadAlgorithm