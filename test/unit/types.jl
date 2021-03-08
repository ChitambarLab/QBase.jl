using Test, QBase

@testset "Testing QBase.jl types" begin
    println("Testing QBase.jl types : ")
    for test in readdir("./test/unit/types/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./unit/types/$test")
            include("./types/$test")
        end
    end
end
