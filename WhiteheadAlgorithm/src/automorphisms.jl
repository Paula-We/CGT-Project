struct Endomorphism 
    A::FreeAlphabet
    images::Vector{Vector{Int}} #Images on the generators (but not on their inverses)
   
    function Endomorphism(A::FreeAlphabet, images::Vector{Vector{Int}}, safe=true::Bool)
        if safe
            @assert length(A)==2*length(images)
            for i in 1:length(images)
                im = images[i]
                for l in im
                    @assert l <= length(A) "$l is too big for the given alphabet"
                end
                images[i]=freerewrite(A, im)
            end
        end
        return new(A, images)
    end
end

function compose(ϕ::Endomorphism,ψ::Endomorphism)
    @assert ϕ.A==ψ.A "The morphisms have to have the same Alphabet."
    A = ϕ.A
    newimages = Vector{Vector{Int}}(undef,length(ϕ.images))
    for i in 1:length(newimages)
        newimages[i]=Vector{Int}()
        for l in ψ.images[i]
            if length(A.letterToGen) < l
                println(l)
                println(A.letterToGen)
                println(length(A))
            end
            (k, isgen) = A.letterToGen[l]
            if isgen #l is a generator itself and not the inverse of a generator
                append!(newimages[i], ϕ.images[k])
            else #l is the inverse of a generator
                append!(newimages[i], inv(ϕ.A, ϕ.images[k]))
            end
        end
    end
    return Endomorphism(ϕ.A, newimages)
end

function Base.show(io::IO, ϕ::Endomorphism)
    l=length(ϕ.images)
    println(io, "Endomorphism of the free group on $l generators")
    for i in 1:length(ϕ.images)
        w=ϕ.images[i]
        im = word(ϕ.A, w) 
        gen = ϕ.A.gen[i]
        println(io, "$gen maps to $im")
    end
end

function isIdentity(ϕ::Endomorphism)
    for i in length(ϕ.images)
        if ϕ.images[i]!= [ϕ.A.gen[i]]
            return false
        end
    end
    return true
end

function evaluate(A::FreeAlphabet, ϕ::Endomorphism, word::Vector{Int})
    image = Vector{Int}()
    for l in word
        (k, isgen) = A.letterToGen[l]
        if isgen #l is a generator itself and not the inverse of a generator
            append!(image, ϕ.images[k])
        else #l is the inverse of a generator
            append!(image, inv(A, ϕ.images[k]))
        end
    end
    return freerewrite(A, image)
end


#This is the Nielsen automorphism which fixes every generator but the ith and sends the ith to a product with the jth
#or with the inverse of the jth (depending on lr and pm)
function NielsenAut(A::FreeAlphabet, i::Int, j::Int, lr::Char, pm::Char)
    @assert i != j "This is only a Nielsen automorphism for $i not equal to $j"
    @assert i <= length(A)/2 && j <=length(A)/2 "$i or $j is too big"
    @assert lr ∈ ['l','r'] "$lr is an invalid argument, only 'l' and 'r' are allowed"
    @assert pm ∈ ['+','-'] "$pm is an invalid argument, only '+' and '-' are allowed"
    images = Vector{Vector{Int}}(undef,floor(Int,length(A)/2))
    for k in 1:length(images)
        if k == i
            if lr == 'l'
                if pm == '+'
                    images[k] = [A.gen[j], A.gen[k]]
                else
                    images[k] = [inv(A, A.gen[j]), A.gen[k]]
                end
            else
                if pm == '+'
                    images[k] = [A.gen[k], A.gen[j]]
                else
                    images[k] = [A.gen[k], inv(A, A.gen[j])]
                end
            end
        else
            images[k]=[A.gen[k]]
        end
    end
    return Endomorphism(A, images, false)
end

#Whiteheadautomorphism that fixes the ith generator a, and sends every other generator x to one of the following options
#1: x^-1
#2: xa
#3: a^⁻1x
#4: a^-1xa
#5: xa^-1
#6: ax
#7: axa^-1
function WhiteheadAut(A::FreeAlphabet, i::Int, options::Vector{Int})
    @assert length(A)==(length(options)+1)*2 "Optionsize isn't matching to Alphabetsize"
    @assert i <= length(A)/2
    images = Vector{Vector{Int}}(undef, length(options)+1)
    images[i]=[A.gen[i]]
    for k in 1:length(options)
        @assert options[k] ∈ 1:7
        skip = k >= i ? 1 : 0
        k2 = k + skip
        if options[k]==1
            images[k2]=[inv(A, A.gen[k2])]
        elseif options[k]==2
            images[k2]=[A.gen[k2], A.gen[i]]
        elseif options[k]==3
            images[k2]=[inv(A, A.gen[i]), A.gen[k2]]
        elseif options[k]==4
            images[k2]=[inv(A, A.gen[i]), A.gen[k2], A.gen[i]]
        elseif options[k]==5
            images[k2]=[A.gen[k2],inv(A, A.gen[i])]
        elseif options[k]==6
            images[k2]=[A.gen[i], A.gen[k2]]
        else
            images[k2]=[A.gen[i], A.gen[k2], inv(A, A.gen[i])]
        end
    end
    return Endomorphism(A, images, false)
end

#WhiteheadAut(A, i, options) is inverse to WhiteheadAut(A, i, inverseoptions(options))
function inverseoptions(options::Vector{Int})
    newoptions = Vector{Int}(undef, length(options))
    v = [1,5,6,7,2,3,4]
    for k in 1:length(options)
        newoptions[k]=v[options[k]]
    end
    return newoptions
end

#Automorphism that permutes the ith and the jth generator
function TranspositionAut(A::FreeAlphabet, i::Int, j::Int)
    @assert i<=length(A)/2 && j<=length(A)/2 "One of the indices is too big for this alphabet"
    images = Vector{Vector{Int}}(undef,floor(Int,length(A)/2))
    for k in 1:length(images)
        if k == i
            images[k]= [A.gen[j]]
        elseif k == j
            images[k]=[A.gen[i]]
        else
            images[k]=[A.gen[k]]
        end
    end
    return Endomorphism(A, images, false)
end

#Automorphism that inverts the ith generator
function InvertGenAut(A::FreeAlphabet, i::Int)
    @assert i<=length(A)/2 "$i is too big for this alphabet"
    n = floor(Int,length(A)/2)
    images = [[A.gen[k]] for k in 1:n]
    images[i] = [inv(A, A.gen[i])]
    return Endomorphism(A, images, false)
end