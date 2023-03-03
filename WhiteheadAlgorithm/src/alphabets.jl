struct Alphabet{T}
    letters :: Vector{T}
	dict :: Dict{T,Int}
	inverse :: Vector{Int}
	function Alphabet{T}(vect::Vector{T}) where T
		@assert !(T <: Integer) "The type cannot be an Integer"
		alphdict = Dict{T, Int}()
		inverse = zeros(length(vect))
		for i in 1:length(vect)
			alphdict[vect[i]]=i
		end
		return new{T}(vect, alphdict, inverse)
	end
end

struct FreeAlphabet{T} 
    letters :: Vector{T}
	dict :: Dict{T,Int}
	inverse :: Vector{Int}
    gen :: Vector{Int} # gen[i] gives us the ith generator of this alphabet
    # letterToGen[i] gives us the index of the generator corresponding to the ith letter
    # and the boolean indicates if this letter is the generator or the inverse of the generator
    letterToGen :: Vector{Tuple{Int,Bool}} 
	function FreeAlphabet(A::Alphabet{T}) where T
		alreadyseen = zeros(Bool, length(A))
        gen = Vector{Int}()
        letterToGen = Vector{Tuple{Int,Bool}}(undef, length(A))
        for i in 1:length(A)           
            if !alreadyseen[i]
                @assert hasinverse(A, i) && inv(A, i) != i #every letter needs another letter as inverse
                alreadyseen[inv(A,i)] = 1
                push!(gen, i)
                letterToGen[i] = (length(gen), true)
                letterToGen[inv(A,i)] = (length(gen), false)
            end
        end
        @assert length(gen)*2 == length(A)
		return new{T}(A.letters, A.dict, A.inverse, gen, letterToGen)
	end
end

Alph = Union{Alphabet, FreeAlphabet}

Base.length(A::Alph) = length(A.letters)
Base.getindex(A::Alph, letter) = A.dict[letter]
Base.getindex(A::Alph, index::Integer) = A.letters[index]

function hasinverse(A::Alph, index::Integer)
    return A.inverse[index] != 0
end

hasinverse(A::Alph, letter::T) where T = hasinverse(A, A[letter])

function Base.inv(A::Alph, index::Integer)
    hasinverse(A, index) ||
        throw(ArgumentError("Non-invertible letter: $(A[index])"))

    return A.inverse[index]
end

function setinverse!(A::Alph, x::Integer, X::Integer)
    @assert !hasinverse(A, x) "Letter $(A[x]) already has inverse: $(inv(A, x))"
    @assert !hasinverse(A, X) "Letter $(A[X]) already has inverse: $(inv(A, X))"

    A.inverse[x] = X
    A.inverse[X] = x

    return A
end

Base.iterate(A::Alph) = iterate(A.letters)
Base.iterate(A::Alph, state) = iterate(A.letters, state)

function Base.show(io::IO, A::Alph)
    println(io, "Alphabet of $(eltype(A.letters[1])) with $(length(A)) letters:")
    for letter in A
        print(io, A[letter], ".\t", letter)
        if hasinverse(A, letter)
            print(io, " with inverse ", A[inv(A, A[letter])])
        end
        println(io, "")
    end
end

#takes a word represented by a vector of integers and makes a string out of it
function word(A::Alph, v::Vector{Int})
    str=""
    for i in v
        str = str*string(A[i])
    end
    return str
end

#generates an alphabet for the free group in n<=26 generators
function generateFreeAlphabet(n::Int) 
    @assert n <= 26 "This function only works up to n = 26"
    alph = ['a', 'A', 'b', 'B', 'c', 'C', 'd', 'D', 'e', 'E', 'f', 'F', 'g', 'G', 'h', 'H', 'i', 'I', 'j', 'J', 'k', 'K', 'l', 'L', 'm', 'M', 'n', 'N', 'o', 'O', 'p', 'P', 'q', 'Q', 'r', 'R', 's', 'S', 't', 'T', 'u', 'U', 'v', 'V', 'w', 'W', 'x', 'X', 'y', 'Y', 'z', 'Z']
    A = Alphabet{Char}(alph[1:2*n])
    for i in 1:n
        setinverse!(A,2*i-1,2*i)
    end
    return FreeAlphabet(A)
end

#generates an alphabet for the free group
function generateBigFreeAlphabet(n::Int)
    letters = Vector{String}(undef, 2*n)
    for i in 1:n
        letters[2*i-1]="x$i"
        letters[2*i]="X$i"
    end
    A = Alphabet{String}(letters)
    for i in 1:n
        setinverse!(A,2*i-1,2*i)
    end
    return FreeAlphabet(A)
end

function Base.inv(A::Alph, word::Vector{Int})
    inverse = Vector{Int}()
    if length(word)==0
        return word
    end
    for i in length(word):-1:1
        push!(inverse, inv(A, word[i]))
    end
    return inverse
end

function summedweight(words::Vector{Vector{Int}})
    weights = length.(words)
    return sum(weights)
end
