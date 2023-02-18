struct Alphabet{T}
    letters :: Vector{T}
	dict :: Dict{T,Int}
	inverse :: Vector{Int}
	function Alphabet{T}(vect) where T
		@assert !(T <: Integer) "The type cannot be an Integer"
		alphdict = Dict{T, Int}()
		inverse = zeros(length(vect))
		for i in 1:length(vect)
			alphdict[vect[i]]=i
		end
		return new{T}(vect, alphdict, inverse)
	end
end

Base.length(A::Alphabet) = length(A.letters)
Base.getindex(A::Alphabet, letter) = A.dict[letter]
Base.getindex(A::Alphabet, index::Integer) = A.letters[index]

function hasinverse(A::Alphabet, index::Integer)
    return A.inverse[index] != 0
end

hasinverse(A::Alphabet, letter::T) where T = hasinverse(A, A[letter])

function Base.inv(A::Alphabet, index::Integer)
    hasinverse(A, index) ||
        throw(ArgumentError("Non-invertible letter: $(A[index])"))

    return A.inverse[index]
end
Base.inv(A::Alphabet, letter::T) where T = A[Base.inv(A,A[letter])]

function setinverse!(A::Alphabet, x::Integer, X::Integer)
    @assert !hasinverse(A, x) "Letter $(A[x]) already has inverse: $(inv(A, x))"
    @assert !hasinverse(A, X) "Letter $(A[X]) already has inverse: $(inv(A, X))"

    A.inverse[x] = X
    A.inverse[X] = x

    return A
end

Base.length(A::Alphabet) = length(A.letters)
Base.iterate(A::Alphabet) = iterate(A.letters)
Base.iterate(A::Alphabet, state) = iterate(A.letters, state)
Base.eltype(::Type{Alphabet{T}}) where {T} = T

function Base.show(io::IO, A::Alphabet)
    println(io, "Alphabet of $(eltype(A)) with $(length(A)) letters:")
    for letter in A
        print(io, A[letter], ".\t", letter)
        if hasinverse(A, letter)
            print(io, " with inverse ", inv(A, letter))
        end
        println(io, "")
    end
end