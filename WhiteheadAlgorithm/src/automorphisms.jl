#One should use these functions with an alphabet from generateFreeAlphabet(), otherwise the generators may be ordered incorrectly
#Warning: i and j refer to the ith and jth generator which are letters 2i-1 and 2j-1 in our alphabet while l is really the lth letter in A
#This function evaluates the Nielsen automorphism which fixes every generator but the ith and sends the ith to a product with the jth
#or with the inverse of the jth
function nielsenonletter(A::Alphabet, i::Int, j::Int, lr::Char, pm::Char, l::Int)
    @assert l <= length(A) "$l not in Alphabet"
    if l == 2i-1
        if lr == 'l'
            if pm == '+'
                return [2j-1,l]
            else
                return [2j,l]
            end
        else
            if pm == '+'
                return [l, 2j-1]
            else
                return [l, 2j]
            end
        end
    elseif l == 2i
        if lr == 'l'
            if pm == '+'
                return [l, 2j]
            else
                return [l, 2j-1]
            end
        else
            if pm == '+'
                return [2j, l]
            else
                return [2j-1, l]
            end
        end
    else
        return [l]
    end
end

function nielsen(A::Alphabet, i::Int, j::Int, lr::Char, pm::Char, word::Vector{Int})
    @assert i != j "This is only a Nielsen automorphism for $i not equal to $j"
    @assert i <= length(A)/2 && j <=length(A)/2 "$i or $j is not in the Alphabet"
    @assert lr ∈ ['l','r'] "$lr is an invalid argument, only 'l' and 'r' are allowed"
    @assert pm ∈ ['+','-'] "$pm is an invalid argument, only '+' and '-' are allowed"
    newword = Vector{Int}()
    if length(word)==0
        return word
    end
    for l in word
        append!(newword, nielsenonletter(A,i,j,lr,pm,l))
    end
    return freerewrite(A, newword)
end

struct Endomorphism 
    A::Alphabet
    images::Vector{Vector{Int}} #Images on the generators (but not on their inverses)
   
    function Endomorphism(A, images)
        @assert length(A)==2*length(images)
        for i in 1:length(images)
            im = images[i]
            for l in im
                @assert l <= 2*length(images) "This does not define an endomorphism"
            end
            images[i]=freerewrite(A,im)
        end
        return new(A, images)
    end
end

function compose(ϕ::Endomorphism,ψ::Endomorphism)
    @assert ϕ.A==ψ.A "The morphisms have to have the same Alphabet."
    newimages = Vector{Vector{Int}}(undef,length(ϕ.images))
    for i in 1:length(newimages)
        newimages[i]=Vector{Int}()
        for l in ψ.images[i]
            if !iseven(l) #l is a generator itself and not the inverse of a generator
                append!(newimages[i], ϕ.images[floor(Int, (l+1)/2)])
            else #l is the inverse of a generator
                append!(newimages[i], inv(ϕ.A, ϕ.images[floor(Int,l/2)]))
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
        gen = ϕ.A[2*i-1]
        println(io, "$gen maps to $im")
    end
end

function isIdentity(ϕ::Endomorphism)
    for i in length(ϕ.images)
        if ϕ.images[i]!= [2i-1]
            return false
        end
    end
    return true
end

function evaluate(A::Alphabet, ϕ::Endomorphism, word::Vector{Int})
    image=Vector{Int}()
    for l in word
        if !iseven(l) #l is a generator itself and not the inverse of a generator
            append!(image, ϕ.images[floor(Int, (l+1)/2)])
        else #l is the inverse of a generator
            append!(image, inv(A, ϕ.images[floor(Int,l/2)]))
        end
    end
    return image
end

function NielsenAut(A, i, j, lr, pm)
    @assert i != j "This is only a Nielsen automorphism for $i not equal to $j"
    @assert i <= length(A)/2 && j <=length(A)/2 "$i or $j is not in the Alphabet"
    @assert lr ∈ ['l','r'] "$lr is an invalid argument, only 'l' and 'r' are allowed"
    @assert pm ∈ ['+','-'] "$pm is an invalid argument, only '+' and '-' are allowed"
    images = Vector{Vector{Int}}(undef,floor(Int,length(A)/2))
    for k in 1:length(images)
        if k == i
            if lr == 'l'
                if pm == '+'
                    images[k] = [2*j-1,2*k-1]
                else
                    images[k] = [2*j,2*k-1]
                end
            else
                if pm == '+'
                    images[k] = [2*k-1,2*j-1]
                else
                    images[k] = [2*k-1,2*j]
                end
            end
        else
            images[k]=[2*k-1]
        end
    end
    return Endomorphism(A, images)
end

#Whiteheadautomorphism that fixes the ith generator a, and sends every other generator x to one of the following options
#1: x^-1
#2: xa
#3: a^⁻1x
#4: a^-1xa
#5: xa^-1
#6: ax
#7: axa^-1
function WhiteheadAut(A::Alphabet, i::Int, options::Vector{Int})
    @assert length(A)==length(options)*2 "Optionsize isn't matching to Alphabetsize"
    @assert i <= length(A)/2
    images = Vector{Vector{Int}}(undef, length(options))
    for k in 1:length(options)
        @assert options[k] ∈ 1:7
        if options[k]==1
            images[k]=[2k]
        elseif options[k]==2
            images[k]=[2k-1,2i-1]
        elseif options[k]==3
            images[k]=[2i,2k-1]
        elseif options[k]==4
            images[k]=[2i,2k-1,2i-1]
        elseif options[k]==5
            images[k]=[2k-1,2i]
        elseif options[k]==6
            images[k]=[2i-1,2k-1]
        else
            images[k]=[2i-1,2k-1,2i]
        end
        if k == i
            images[k]=[2k-1]
        end
    end
    return Endomorphism(A, images)
end

#WhiteheadAut(A, i, options) is inverse to WhiteheadAut(A, i inverseoptions(options))
function inverseoptions(options::Vector{Int})
    newoptions = Vector{Int}(undef, length(options))
    for k in 1:length(options)
        if options[k]==1
            newoptions[k]=1
        elseif options[k]==2
            newoptions[k]=5
        elseif options[k]==3
            newoptions[k]=6
        elseif options[k]==4
            newoptions[k]=7
        elseif options[k]==5
            newoptions[k]=2
        elseif options[k]==6
            newoptions[k]=3
        else
            newoptions[k]=4
        end
    end
    return newoptions
end
