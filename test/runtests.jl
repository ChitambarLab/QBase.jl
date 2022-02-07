using Test, SafeTestsets

println("importing QBase.jl")
@time using QBase

println("Testing QBase.jl")
@time @testset "./src/QBase.jl"  begin

println("Testing Math Utilities : ")
@time @testset "Math Utilities" begin
    println("./test/unit/math/matrices.jl")
    @time @safetestset "./test/unit/math/matrices.jl" begin include("unit/math/matrices.jl") end

    println("./test/unit/math/combinatorics.jl")
    @time @safetestset "./test/unit/math/combinatorics.jl" begin include("unit/math/combinatorics.jl") end
end

println("Testing Types : ")
@time @testset "Types" begin
    println("./test/unit/types/probabilities.jl")
    @time @safetestset "./test/unit/types/probabilities.jl" begin include("unit/types/probabilities.jl") end

    println("./test/unit/types/brakets.jl")
    @time @safetestset "./test/unit/types/brakets.jl" begin include("unit/types/brakets.jl") end

    println("./test/unit/types/states.jl")
    @time @safetestset "./test/unit/types/states.jl" begin include("unit/types/states.jl") end

    println("./test/unit/types/unitaries.jl")
    @time @safetestset "./test/unit/types/unitaries.jl" begin include("unit/types/unitaries.jl") end

    println("./test/unit/types/measurements.jl")
    @time @safetestset "./test/unit/types/measurements.jl" begin include("unit/types/measurements.jl") end

    println("./test/unit/kraus_channels.jl")
    @time @safetestset "./test/unit/kraus_channels.jl" begin include("unit/kraus_channels.jl") end

    println("./test/unit/choi_operators.jl")
    @time @safetestset "./test/unit/choi_operators.jl" begin include("unit/choi_operators.jl") end
end

println("Testing Constructors : ")
@time @testset "Constructors" begin
    println("./test/unit/constructors/brakets.jl")
    @time @safetestset "./test/unit/constructors/brakets.jl" begin include("unit/constructors/brakets.jl") end

    println("./test/unit/constructors/states.jl")
    @time @safetestset "./test/unit/constructors/states.jl" begin include("unit/constructors/states.jl") end

    println("./test/unit/constructors/unitaries.jl")
    @time @safetestset "./test/unit/constructors/unitaries.jl" begin include("unit/constructors/unitaries.jl") end

    println("./test/unit/constructors/measurements.jl")
    @time @safetestset "./test/unit/constructors/measurements.jl" begin include("unit/constructors/measurements.jl") end
end

println("Testing Main : ")
@time @testset "Main" begin
    println("./test/unit/evolve.jl")
    @time @safetestset "./test/unit/evolve.jl" begin include("unit/evolve.jl") end

    println("./test/unit/measure.jl")
    @time @safetestset "./test/unit/measure.jl" begin include("unit/measure.jl") end

    println("./test/unit/channels.jl")
    @time @safetestset "./test/unit/channels.jl" begin include("unit/channels.jl") end

    println("./test/unit/information.jl")
    @time @safetestset "./test/unit/information.jl" begin include("unit/information.jl") end
end

end
