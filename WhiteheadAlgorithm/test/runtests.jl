using Test
using WhiteheadAlgorithm

WHA = WhiteheadAlgorithm
@testset "alphabets" begin
    A = WHA.Alphabet{Char}(['x','y','z'])
    WHA.setinverse!(A,1,3)
    @test WHA.inv(A, 3) == 1
    @test WHA.length(A) == 3
    @test A[1] == 'x'
    @test A['x'] == 1
    @test WHA.hasinverse(A, 1) == true
    @test WHA.hasinverse(A, 2) == false
    @test_throws ArgumentError WHA.inv(A,2)
    @test_throws AssertionError WHA.setinverse!(A, 2, 3)
end

@testset "free rewrite" begin
    A = WHA.Alphabet{Char}(['x','X','y','Y'])
    WHA.setinverse!(A,1,2)
    WHA.setinverse!(A,3,4)
    @test WHA.freerewrite(A, Vector{Int}()) == []
    @test WHA.freerewrite(A, [1,2]) == []
    @test WHA.freerewrite(A, [1,4,3,2,2]) == [2]
    @test WHA.freerewrite(A, [1,3,2,4]) == [1,3,2,4]
    @test_throws AssertionError WHA.freerewrite(A, [5])
end

@testset "Nielsen automorphisms" begin
    A = WHA.generateFreeAlphabet(3)
    @test WHA.evaluate(WHA.NielsenAut(A,1,2,'l','-'), [1]) == [4,1]
    @test WHA.evaluate(WHA.NielsenAut(A,1,2,'r','+'), [1]) == [1,3]
    @test WHA.evaluate(WHA.NielsenAut(A,1,2,'r','-'), [1]) == [1,4]
    @test_throws AssertionError WHA.NielsenAut(A,1,1,'l','+')
    N = WHA.NielsenAut(A,1,2,'l','+')
    @test (w -> WHA.evaluate(N, w)).([[1],[3],[5]])==[[3,1],[3],[5]] 
end

@testset "WhiteheadAut" begin
    A = WHA.generateFreeAlphabet(8)
    options = [2,3,4,5,6,7,1]
    WAut = WHA.WhiteheadAut(A, 1, options)
    invopt = WHA.inverseoptions(options)
    WAut2 = WHA.WhiteheadAut(A, 1, invopt)
    @test WHA.isIdentity(WHA.compose(WAut, WAut2))==true
    @test (w -> WHA.evaluate(WAut, w)).([[1],[3],[5],[7],[9],[11],[13],[15]])==[[1],[3,1],[2,5],[2,7,1],[9,2],[1,11],[1,13,2],[15]]
end

@testset "Part of Basis" begin
    A = WHA.generateFreeAlphabet(3)
    w = [3,1,4,6,2,3]
    (b, f) = WHA.isPrimitive(A, w,false) #here we should need both additional transformations
    @test b == true
    @test WHA.evaluate(f, [1]) == w
    (b, f) = WHA.isPrimitive(A, [2,5,4,6,1,3,6,3,1],false) #here we also need Whitehead automorphisms
    @test b == true
    (b, f) = WHA.isPrimitive(A, [2,5,4,6,1,3,1,3,4,2,6,3,1,5],false)
    @test b == false
    (b, f) = WHA.isPartOfBasis(A, [[1,1]],false)
    @test b == false
    (b, f) = WHA.isPartOfBasis(A, [[1,3,2,4]],false)
    @test b == false   
    (b, f) = WHA.isPartOfBasis(A, [[1,2]],false) #Trivial element is not part of basis
    @test b == false
    (b, f) = WHA.isPartOfBasis(A, [[1],[2,4,1],[1]],false) #All elements in a basis have to be pairwise distinct
    @test b == false    
    (b, f) = WHA.isPartOfBasis(A, [[2,5,4,6,1,3,6,3,1],[2,4,1],[2,5,4,1]],false)
    @test b == true
    @test WHA.evaluate(f, [1]) == [2,5,4,6,1,3,6,3,1]
    @test WHA.evaluate(f, [3]) == [2,4,1]
    @test WHA.evaluate(f, [5]) == [2,5,4,1]
end