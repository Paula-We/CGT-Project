module WhiteheadAlgorithm
WHA = WhiteheadAlgorithm

include("alphabets.jl")
include("freerewrite.jl")
include("automorphisms.jl")

function isPrimitive(A::Alphabet, w0::Vector{Int}, pr=true::Bool)
    if length(freerewrite(A,w0))==0
        println("The trivial element is not primitive")
        return (false, f)
    end
    (v,f)=minimize(A,[w0],pr)
    w=v[1]
    if length(w)==1 #w0 is primitive
        (v,f)=permuteGensandInv(A,v,f,pr)
        w=v[1]
        pr&&println(word(A, w0)*" is a primitive word. Here is an automorphism mapping the first generator to it:")
        pr&&print(f)
        return (true,f)
    else
        pr&&println(word(A, w0)*" is not a primitive word, the minimal element in its equivalence class is "*word(A, w))
        return (false, f)
    end
end

function isPartOfBasis(A::Alphabet, b0::Vector{Vector{Int}}, pr=true::Bool) 
    (b,f)=minimize(A,b0,pr)
    if summedweight(b)==length(b)
        (b,f)=permuteGensandInv(A,b,f,pr)
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
function minimize(A::Alphabet, b0::Vector{Vector{Int}}, pr=true::Bool)
    b = (w -> freerewrite(A,w)).(b0)
    pr&&println((w -> word(A,w)).(b0))
    pr&&println("")
    n = floor(Int,length(A)/2)
    f = Endomorphism(A, [[2k-1] for k in 1:n]) #f is initialized as Identity
    autfound = true
    for w in b
        if length(w)==0
            println("The trivial element can not be part of a basis")
            return (false, f)
        end
    end
    @assert length(b0) <= n "Too many elements to be a basis of the free group in $n generators"
    #repeat till no length shortening automorphism is found or every word in the vector has length 1
    while autfound && summedweight(b)!=length(b) 
        autfound = false
        #iterate over all Nielsen automs except for x |-> x^-1 because they don't shorten any word
        for (i,j,lr,pm) in Iterators.product(1:n,1:n,['l','r'],['+','-']) 
            if i != j
                b2 = (w -> nielsen(A,i,j,lr,pm,w)).(b)
                if summedweight(b2) < summedweight(b) #the automorphism shortened b
                    b = b2
                    pr&&println(NielsenAut(A,i,j,lr,pm))
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
        if autfound == true
            continue
        end
        pr&&println("Now we need to look at all Whitehead automorphisms...")
        for i in 1:n
            #iterate over all Whitehead automorphisms except for the permutations because they don't change the length
            for options in Iterators.product([collect(1:4) for k in 1:(n-1)]...) 
                options = collect(options)
                WAut = WhiteheadAut(A, i, options)
                b2 = (w -> evaluate(A, WAut, w)).(b)
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
                b3 = (w -> evaluate(A, WAut2, w)).(b)
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

#This function takes a vector of words where each word is either a generator or an inverse of a generator
#and finds an automorphism taking the words to the first generators
#We also compose this automorphism with the the one we already had before
function permuteGensandInv(A::Alphabet, b::Vector{Vector{Int}}, f::Endomorphism, pr=true::Bool)
    pr&&println("")
    #First we map our word vector to another one where each word is a generator 
    #(so no inverses of generators anymore)
    for i in 1:length(b)
        if iseven(b[i][1])
            Aut = InvertGenAut(A, floor(Int, b[i][1]/2))
            b = (w -> evaluate(A, Aut, w)).(b)
            f = compose(f, Aut)
            pr&&println(Aut)
            pr&&println((w -> word(A,w)).(b))
            pr&&println("")
        end
    end
    #Now we need to get the generators in the right ordered
    #To achieve that we use transpositions
    for j in 1:length(b)
        i = floor(Int, (b[j][1]+1)/2) #the b[j]th-letter of our alphabet corresponds to the ith generator
        if i != j
            T = TranspositionAut(A, i, j)
            b = (w -> evaluate(A, T, w)).(b)
            f = compose(f, T)
            pr&&println(T)
            pr&&println((w -> word(A,w)).(b))
            pr&&println("")
        end
    end
    return (b,f)
end

#[3,1,4,6,2,3] is primitive (only using Nielsen automorphisms) in F_3 (here we need both additional transformations at the end)
#[1,1] is not primitive in F_2
#[2,5,4,6,1,3,1,3,4,2,6,3,1,5] not primitive in F_3
#[2,5,4,6,1,3,6,3,1] is primitive also using Whitehead automorphisms in F_3
#[[2,5,4,6,1,3,6,3,1],[2,4,1],[2,5,4,1]] is a basis
A = generateFreeAlphabet(5)
isPartOfBasis(A, [[1,3,2,4,4,2,3,1]])
#isPrimitive(A, [3,1,4,6,2,3])

end # module WhiteheadAlgorithm