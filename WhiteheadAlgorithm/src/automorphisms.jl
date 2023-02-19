#One should use these functions with an alphabet from generateFreeAlphabet(), otherwise the generators may be ordered incorrectly
#Warning: i and j refer to the ith and jth generator which are letters 2i-1 and 2j-1 in our alphabet while l is really the lth letter in A
#This function evaluates the Nielsen automorphism which fixes every generator but the ith and sends the ith to a product with the jth
#or with the inverse of the jth
function nielsenonletter(A::Alphabet, i::Int, j::Int, lr::Char, pm::Char, l::Int)
    @assert i != j "This is only a Nielsen automorphism for $i not equal to $j"
    @assert i <= length(A)/2 && j <=length(A)/2 "$i or $j is not in the Alphabet"
    @assert lr ∈ ['l','r'] "$lr is an invalid argument, only 'l' and 'r' are allowed"
    @assert pm ∈ ['+','-'] "$pm is an invalid argument, only '+' and '-' are allowed"
    @assert l <= length(A)
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
                append!(newimages[i], inv(A, ϕ.images[floor(Int,l/2)]))
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


