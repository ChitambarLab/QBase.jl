using Test, SafeTestsets

println("importing QBase.jl")
@time using QBase

function _test_runner()
    @time @safetestset "./test/unit/constructors/states.jl" begin include("unit/constructors/states.jl") end


    @testset "QBase.jl" begin
        @testset "unit tests:" begin
            println("running unit tests.")
            for test in readdir("./test/unit/")
                # run only julia files in test directory
                if occursin(r"^.*\.jl$", test)
                    println("./unit/$test")
                    @time include("./unit/$test")
                end
            end

            println("Testing math utilities : ")
            for test in readdir("./test/unit/math/")
                # run only julia files in test directory
                if occursin(r"^.*\.jl$", test)
                    println("./unit/math/$test")
                    @time include("unit/math/$test")
                end
            end

            println("Testing QBase.jl types : ")
            # for test in readdir("./test/unit/types/")
            #     # run only julia files in test directory
            #     if occursin(r"^.*\.jl$", test)
            #         println("./unit/types/$test")
            #         @time include("unit/types/$test")
            #     end
            # end
            @time @safetestset "./test/unit/types/brakets.jl" begin include("unit/types/brakets.jl") end
            @time @safetestset "./test/unit/types/states.jl" begin include("unit/types/states.jl") end
            @time @safetestset "./test/unit/types/probabilities.jl" begin include("unit/types/probabilities.jl") end

            println("Testing constructors : ")
            @time @safetestset "./test/unit/constructors/states.jl" begin include("unit/constructors/states.jl") end
            @time @safetestset "./test/unit/constructors/brakets.jl" begin include("unit/constructors/brakets.jl") end


        end
    end


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
