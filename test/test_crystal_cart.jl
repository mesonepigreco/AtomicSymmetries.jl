using AtomicSymmetries
using LinearAlgebra
using Test


function test_crystal_to_cart(; verbose=false)
    # Let us define a wired primitive cell
    cell = [1.0 0.2 0.5; 0.2 1.0 0.2; 0.2 0.2 1.0]

    q_points_fract = [0.0 0.5 0.25 0.25;
                      0.0 0.5 0.5 0.25;
                      0.0 0.5 0.5 0.5]

    reciprocal_lattice = zeros(Float64, 3, 3)
    get_reciprocal_lattice!(reciprocal_lattice, cell)

    if verbose
        @show reciprocal_lattice
    end

    q_points_cartesian = zeros(Float64, size(q_points_fract)...)
    q_points_cartesian_test = zeros(Float64, size(q_points_fract)...)
    q_points_cryst_test = zeros(Float64, size(q_points_fract)...)

    mul!(q_points_cartesian, reciprocal_lattice, q_points_fract, 1/(2π), 0.0)
    mul!(q_points_cryst_test, cell', q_points_cartesian, 1.0, 0.0)

    # Check if the two conversion went well
    @test isapprox(q_points_cryst_test, q_points_fract; atol = 1e-8, rtol =1e-6)

    # Test with cryst_cart_conv!
    cryst_cart_conv!(q_points_cartesian_test, q_points_fract, cell, reciprocal_lattice,
                     true; q_space=true)

    @test isapprox(q_points_cartesian, q_points_cartesian_test; atol=1e-12, rtol=1e-8)

    cryst_cart_conv!(q_points_cryst_test, q_points_cartesian, cell, reciprocal_lattice,
                     false; q_space=true)

    @test isapprox(q_points_fract, q_points_cryst_test; atol=1e-12, rtol=1e-8)



    # Let us define atomic positions with the following crystal coordinates
    cryst_coords = [0.0 0.8; 0.5 0.2; 0.0 -0.8]

    # Now we convert the crystal coordinates into cartesian ones
    cart_coords = zeros(Float64, 3, 2)
    new_cart_coords = zeros(Float64, 3, 2)
    new_cryst_coords = zeros(Float64, 3, 2)

    get_cartesian_coords!(cart_coords, cryst_coords, cell)
    get_crystal_coords!(new_cryst_coords, cart_coords, cell)

    if verbose
        @show cryst_coords'
        @show cart_coords'
        @show new_cryst_coords'
    end

    for i in 1:2
        for j in 1:3
            @test new_cryst_coords[j, i] ≈ cryst_coords[j, i] rtol = 1e-8 atol = 1e-12
        end
    end

    # Try again with cryst_cart_conv!
    cryst_cart_conv!(new_cart_coords, cryst_coords, cell, reciprocal_lattice,
                     true; q_space = false)

    @test isapprox(cart_coords, new_cart_coords, rtol=1e-8, atol=1e-12)

    cryst_cart_conv!(new_cryst_coords, cart_coords, cell, reciprocal_lattice,
                     false; q_space = false)

    @test isapprox(cryst_coords, new_cryst_coords, rtol=1e-8, atol=1e-12)


end

if abspath(PROGRAM_FILE) == @__FILE__
    test_crystal_to_cart(; verbose=true)
end

