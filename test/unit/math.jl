using Test, QBase

@testset "Testing math utilites" begin
    println("Testing math utilities : ")
    for test in readdir("./test/unit/math/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./unit/math/$test")
            @time include("./math/$test")
        end
    end
end
