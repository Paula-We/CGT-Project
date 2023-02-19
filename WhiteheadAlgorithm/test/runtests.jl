using Test
using WhiteheadAlgorithm

WHA = WhiteheadAlgorithm
@testset "alphabets" begin
    A = WHA.Alphabet{Char}(['x','y','z'])
    WHA.setinverse!(A,1,3)
    @test WHA.inv(A, 'x') == 'z'
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
    @test WHA.nielsen(A,1,2,'l','+',[1]) == [3,1]
    @test WHA.nielsen(A,1,2,'l','-',[1]) == [4,1]
    @test WHA.nielsen(A,1,2,'r','+',[1]) == [1,3]
    @test WHA.nielsen(A,1,2,'r','-',[1]) == [1,4]
    @test WHA.nielsen(A,1,2,'l','+',[2]) == [2,4]
    @test WHA.nielsen(A,1,2,'l','+',[3]) == [3]
    @test_throws AssertionError WHA.nielsen(A,1,1,'l','+',[1])
    @test WHA.nielsen(A,1,2,'l','+',[4,1]) == [1]
end

@testset "WhiteheadAut" begin
    A = WHA.generateFreeAlphabet(8)
    options = [1,1,2,3,4,5,6,7]
    WAut = WHA.WhiteheadAut(A, 1, options)
    invopt = WHA.inverseoptions(options)
    WAut2 = WHA.WhiteheadAut(A, 1, invopt)
    @test WHA.isIdentity(WHA.compose(WAut, WAut2))==true
    @test WHA.evaluate(A, WAut, [1]) == [1]
    @test WHA.evaluate(A, WAut, [2]) == [2]
    @test WHA.evaluate(A, WAut, [3]) == [4]
    @test WHA.evaluate(A, WAut, [4]) == [3]
    @test WHA.evaluate(A, WAut, [5]) == [5,1]
    @test WHA.evaluate(A, WAut, [6]) == [2,6]
end