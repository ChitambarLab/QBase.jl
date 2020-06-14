using Test, LinearAlgebra

@testset "./src/QBase/Unitaries.jl" begin

using QBase

@testset "QBase.Unitaries.is_unitary()" begin
    valid_unitaries = [
        [0 1;1 0],
        [0 -im 0;im 0 0;0 0 1],
        [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0],
    ]

    num_valid_cases = length(valid_unitaries)

    @testset "valid unitaries $i/$num_valid_cases" for i in 1:num_valid_cases
        @test QBase.is_unitary(valid_unitaries[i])
    end

    invalid_unitaries = [
        [0.5 0;0 0.5],
        [0 im 0;-im 0 0;0 0 0],
    ]

    num_invalid_cases = length(invalid_unitaries)

    @testset "invalid unitaries $i/$num_invalid_cases" for i in 1:num_invalid_cases
        @test !QBase.is_unitary(invalid_unitaries[i])
    end

    @test_throws DomainError QBase.is_unitary([1 0 0;0 1 0])
end

@testset "QBase.Unitaries.Unitary()" begin
    valid_unitaries = [
        [1 0 0;0 1 0;0 0 1],
        [0 1;1 0],
        [1 0 0 0;0 -1 0 0;0 0 0 1;0 0 1 0],
    ]

    @testset "valid unitary: $id/$(length(valid_unitaries))" for id in 1:length(valid_unitaries)
        U = QBase.Unitary(valid_unitaries[id])
        @test isa(U, QBase.Unitary)
        @test U == valid_unitaries[id]
    end

    invalid_unitaries = [
        [0 1;1 0;0 0],
        [1 0;0 2],
        [0.5 0.5;-0.5 0.5]
    ]

    @testset "invalid unitary: $id/$(length(invalid_unitaries))" for id in 1:length(invalid_unitaries)
        U = invalid_unitaries[id]
        @test_throws DomainError QBase.Unitary(U)
    end
end

@testset "QBase.Unitaries.QubitUnitary()" begin
    valid_qubit_unitaries = [
        [0 1;1 0],
        [0 -im;im 0]
    ]

    num_valid = length(valid_qubit_unitaries)
    @testset "valid cases: $i/$num_valid" for i in 1:num_valid
        U = QBase.QubitUnitary(valid_qubit_unitaries[i])
        @test U isa QBase.QubitUnitary
        @test U == valid_qubit_unitaries[i]
    end

    invalid_qubit_unitaries = [
        [1 0 0;0 1 0;0 0 1],
        [1 0;0 0]
    ]

    num_invalid = length(valid_qubit_unitaries)
    @testset "invalid cases: $i/$num_invalid" for i in 1:num_invalid
        @test_throws DomainError QBase.QubitUnitary(invalid_qubit_unitaries[i])
    end
end

@testset "pauli constants" begin
    @test QBase.σx == [0 1; 1 0]
    @test QBase.σx isa QBase.QubitUnitary

    @test QBase.σy == [0 -im; im 0]
    @test QBase.σy isa QBase.QubitUnitary

    @test QBase.σz == [1 0; 0 -1]
    @test QBase.σz isa QBase.QubitUnitary
end

@testset "QBase.Unitaries.paulis" begin
    @test length(QBase.paulis) == 3
    @test QBase.paulis isa Vector{QBase.QubitUnitary}
    @test QBase.paulis[1] == [0 1;1 0]
    @test QBase.paulis[2] == [0 -im;im 0]
    @test QBase.paulis[3] == [1 0;0 -1]
end

@testset "qubit_rotation()" begin

    ρx = 0.5*[1 1; 1 1]
    ρz = [1 0; 0 0]

    @testset "rotation of θ = $θ about x-axis" for θ in -π:0.1π:π
        Rx = QBase.qubit_rotation(θ, axis="x")
        @test isapprox(Rx*ρz*Rx',  [cos(θ/2)^2 im*sin(θ/2)cos(θ/2); -im*sin(θ/2)cos(θ/2) sin(θ/2)^2], atol=1e-7)
    end

    @testset "rotation of θ = $θ about y-axis" for θ in -π:0.1π:π
        Ry = QBase.qubit_rotation(θ, axis="y")
        @test isapprox(Ry*ρz*Ry',  [cos(θ/2)^2 sin(θ/2)cos(θ/2); sin(θ/2)cos(θ/2) sin(θ/2)^2], atol=1e-7)
    end

    @testset "rotation of θ = $θ about z-axis" for θ in -π:0.1π:π
        Rz = QBase.qubit_rotation(θ, axis="z")
        @test isapprox(Rz*ρz*Rz', ρz)

        @test isapprox(Rz*ρx*Rz', [0.5 0.5*(cos(θ)-im*sin(θ)); 0.5*(cos(θ)+im*sin(θ)) 0.5], atol=1e-7)
    end
end

end
