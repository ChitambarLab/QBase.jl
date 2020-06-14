using Test

@testset "running tests for module QBase.jl" begin

    println("running QBase.jl unit tests:")
    for test in readdir("./test/QBase/unit/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./QBase/unit/$test")
            include("./unit/$test")
        end
    end
end
