using Test

@testset "running all tests for module QMath.jl" begin

    println("running QMath.jl unit tests:")
    for test in readdir("./test/unit/QMath/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./unit/QMath/$test")
            include("./QMath/$test")
        end
    end
end
