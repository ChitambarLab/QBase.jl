using Test, SafeTestsets

println("importing QBase.jl")
@time using QBase

@time begin

function _test_runner()

    # @testset "QBase.jl" begin
    #     @testset "unit tests:" begin
    #         println("running unit tests.")
    #         for test in readdir("./test/unit/")
    #             # run only julia files in test directory
    #             if occursin(r"^.*\.jl$", test)
    #                 println("./unit/$test")
    #                 @time include("./unit/$test")
    #             end
    #         end
    #     end
    # end

    println("Testing math utilities : ")
    @time @safetestset "./test/unit/math/matrices.jl" begin include("unit/math/matrices.jl") end
    @time @safetestset "./test/unit/math/combinatorics.jl" begin include("unit/math/combinatorics.jl") end

    println("Testing types : ")
    @time @safetestset "./test/unit/types/probabilities.jl" begin include("unit/types/probabilities.jl") end
    @time @safetestset "./test/unit/types/brakets.jl" begin include("unit/types/brakets.jl") end
    @time @safetestset "./test/unit/types/states.jl" begin include("unit/types/states.jl") end
    @time @safetestset "./test/unit/types/unitaries.jl" begin include("unit/types/unitaries.jl") end
    @time @safetestset "./test/unit/types/measurements.jl" begin include("unit/types/measurements.jl") end

    println("Testing constructors : ")
    @time @safetestset "./test/unit/constructors/brakets.jl" begin include("unit/constructors/brakets.jl") end
    @time @safetestset "./test/unit/constructors/states.jl" begin include("unit/constructors/states.jl") end
    @time @safetestset "./test/unit/constructors/unitaries.jl" begin include("unit/constructors/unitaries.jl") end
    @time @safetestset "./test/unit/constructors/brakets.jl" begin include("unit/constructors/measurements.jl") end

    println("Testing high-level : ")
    @time @safetestset "./test/unit/channels.jl" begin include("unit/channels.jl") end
    @time @safetestset "./test/unit/information.jl" begin include("unit/information.jl") end

end

# Pkg.test("QBase") runs from ./test directory. Development tests from root.
dir = pwd()
if occursin(r".*test$", dir)
    cd(_test_runner, "../")
elseif occursin(r".*QBase", dir)
    _test_runner()
else
    error("runtests.jl is currently running from the $(pwd()) directory with contents $(readdir()). runtests.jl must be run from the ./QBase.jl or ./QBase.jl/test directories.")
end

println("elapsed time :")
end
