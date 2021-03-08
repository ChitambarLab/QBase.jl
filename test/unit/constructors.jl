using Test, QBase

@testset "Testing QBase.jl constructors" begin
    println("Testing QBase.jl constructors : ")
    for test in readdir("./test/unit/constructors/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./unit/constructors/$test")
            include("./constructors/$test")
        end
    end
end
