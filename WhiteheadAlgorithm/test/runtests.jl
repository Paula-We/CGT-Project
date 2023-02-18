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
    @test WHA.freerewrite(A, []) == []
    @test WHA.freerewrite(A, [1,2]) == []
    @test WHA.freerewrite(A, [1,4,3,2,2]) == [2]
    @test WHA.freerewrite(A, [1,3,2,4]) == [1,3,2,4]
    @test_throws AssertionError WHA.freerewrite(A, [5])
end