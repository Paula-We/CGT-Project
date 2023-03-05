#This function takes a word over an alphabet and reduces it freely
function freerewrite(A::Alph, word::Vector{Int})
    redword = Vector{Int}()
    if length(word)==0
        return []
    end
    for i in word
        @assert i <= length(A)
        if length(redword)==0
            push!(redword,i)
        elseif hasinverse(A, last(redword)) && inv(A, last(redword)) == i
            resize!(redword, length(redword)-1)
        else
            push!(redword,i) 
        end
    end
    return redword
end