module WhiteheadAlgorithm
WHA = WhiteheadAlgorithm

include("alphabets.jl")
include("freerewrite.jl")
include("automorphisms.jl")

function isPrimitive(A::Alphabet, w0::Vector{Int}, pr=true::Bool)
    w = freerewrite(A,w0)
    pr&&println(word(A,w))
    pr&&println("")
    n = floor(Int,length(A)/2)
    f = Endomorphism(A, [[2k-1] for k in 1:n]) #f is initialized as Identity
    autfound = true
    if length(w)==0
        println("The trivial element is not primitive")
        return (false, f)
    end
    while autfound && length(w)!=1
        autfound = false
        for (i,j,lr,pm) in Iterators.product(1:n,1:n,['l','r'],['+','-']) #iterate over all Nielsen automs except for x |-> x^-1
            if i != j
                w2 = nielsen(A, i, j, lr, pm, w)
                if length(w2) < length(w)
                    w = w2
                    pr&&println(NielsenAut(A,i,j,lr,pm))
                    pr&&println(word(A,w))
                    pr&&println("")                    
                    pm2 = (pm == '+') ? '-' : '+'
                    f2 = NielsenAut(A, i, j, lr, pm2) #f2 is the inverse of the Nielsen autom used above to shorten w
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
            for options in Iterators.product([collect(1:4) for k in 1:n]...) #iterate over all Whitehead autom except for the permutations
                options = collect(options)
                WAut = WhiteheadAut(A, i, options)
                w2 = evaluate(A, WAut, w)
                invopt = inverseoptions(options)
                WAut2 = WhiteheadAut(A, i, invopt) #WAut2 is the inverse of WAut             
                if length(w2) < length(w)
                    w = w2
                    pr&&println(WAut)
                    pr&&println(word(A,w))
                    f = compose(f,WAut2)
                    autfound = true
                    break
                end
                w3 = evaluate(A, WAut2, w)
                if length(w3) < length(w)
                    w = w3
                    pr&&println(WAut2)
                    pr&&println(word(A,w))
                    f = compose(f,WAut)
                    autfound = true
                    break
                end
            end
        end
    end
    if length(w)==1 #w0 is primitive, now we need two additional automorphisms s.t. the first generator maps to w0
        if iseven(w[1]) #apply automorphism mapping w[1] to w[1]^-1
            images = [[2k-1] for k in 1:n]
            images[floor(Int,w[1]/2)] = [w[1]]
            Autom = Endomorphism(A, images)
            w = evaluate(A, Autom, w)
            f = compose(f, Autom)
        end
        if w[1] != 1 #apply transposition swapping w[1] and 1
            T = TranspositionAut(A, 1, floor(Int, (w[1]+1)/2))
            w = evaluate(A,T,w)
            f = compose(f, T)
        end
        @assert w == [1] "Something went wrong...here"
        @assert f.images[1] == freerewrite(A,w0) "Something went wrong..."
        pr&&println(word(A, w0)*" is a primitive word. Here is an automorphism mapping the first generator to it:")
        pr&&print(f)
        return (true,f)
    else
        pr&&println(word(A, w0)*" is not a primitive word, the minimal element in its equivalence class is "*word(A, w))
        return (false, f)
    end
end

function isPartOfBasis(A::Alphabet, b0::Vector{Vector{Int}}, pr=true::Bool) #We do essentially the same as in isPrimitive, but we have to vectorize some functions
    b = (w -> freerewrite(A,w)).(b0)
    pr&&println((w -> word(A,w)).(b0))
    pr&&println("")
    n = floor(Int,length(A)/2)
    f = Endomorphism(A, [[2k-1] for k in 1:n]) #f is initialized as Identity
    autfound = true
    for w in b
        if length(w)==0
            println("The trivial elemnt can not be part of a basis")
            return (false, f)
        end
    end
    @assert length(b0) <= n "Too many elements to be a basis of the free group in $n generators"
    while autfound && summedweight(b)!=length(b)
        autfound = false
        for (i,j,lr,pm) in Iterators.product(1:n,1:n,['l','r'],['+','-']) #iterate over all Nielsen automs except for x |-> x^-1
            if i != j
                b2 = (w -> nielsen(A,i,j,lr,pm,w)).(b)
                if summedweight(b2) < summedweight(b)
                    b = b2
                    pr&&println(NielsenAut(A,i,j,lr,pm))
                    pr&&println((w -> word(A,w)).(b))
                    pr&&println("")
                    pm2 = (pm == '+') ? '-' : '+'                    
                    f2 = NielsenAut(A, i, j, lr, pm2) #f2 is the inverse of the Nielsen autom used above to shorten b
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
            for options in Iterators.product([collect(1:4) for k in 1:n]...) #iterate over all Whitehead autom except for the permutations
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
    if summedweight(b)==length(b)
        for i in 1:length(b)
            if iseven(b[i][1])
                Aut = InvertGenAut(A, floor(Int, b[i][1]/2))
                b = (w -> evaluate(A, Aut, w)).(b)
                f = compose(f, Aut)
            end
        end
        pr&&println((w -> word(A,w)).(b))
        for j in 1:length(b)
            i = floor(Int, b[j][1]+1/2)
            if i != j
                T = TranspositionAut(A, i, j)
                b = (w -> evaluate(A, T, w)).(b)
                f = compose(f, Aut)
            end
        end
        pr&&println((w -> word(A,w)).(b))
        pr&&println("These elements can be part of a basis. Here's an automorphism mapping some generators to them:")
        pr&&print(f)
        return (true, f)
    else
        pr&&println("These elements can not be part of a basis.")
    end
end
    

#[3,1,4,6,2,3] is primitive (only using Nielsen automorphisms) in F_3 (here we need both additional transformations at the end)
#[1,1] is not primitive in F_2
#[2,5,4,6,1,3,1,3,4,2,6,3,1,5] not primitive in F_3
#[2,5,4,6,1,3,6,3,1] is primitive also using Whitehead automorphisms in F_3
#[[2,5,4,6,1,3,6,3,1],[2,4,1],[2,5,4,1]] is a basis
A = generateFreeAlphabet(3)
isPartOfBasis(A, [[2,5,4,6,1,3,6,3,1],[2,4,1],[2,5,4,1]])
#isPrimitive(A, [2,5,4,6,1,3,6,3,1])

end # module WhiteheadAlgorithm