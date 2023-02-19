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