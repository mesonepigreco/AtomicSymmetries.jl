using AtomicSymmetries
using Test


function test_crystal_to_cart()
    # Let us define a wired primitive cell
    cell = [1.0 0.2 0.5; 0.2 1.0 0.2; 0.2 0.2 1.0]

    # Let us define atomic positions with the following crystal coordinates
    cryst_coords = [0.0 0.8; 0.5 0.2; 0.0 -0.8]

    # Now we convert the crystal coordinates into cartesian ones
    cart_coords = zeros(Float64, 3, 2)
    new_cryst_coords = zeros(Float64, 3, 2)

    get_cartesian_coords!(cart_coords, cryst_coords, cell)
    get_crystal_coords!(new_cryst_coords, cart_coords, cell)

    @show cryst_coords'
    @show cart_coords'
    @show new_cryst_coords'

    for i in 1:2
        for j in 1:3
            @test new_cryst_coords[j, i] ≈ cryst_coords[j, i] rtol = 1e-8 atol = 1e-12
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_crystal_to_cart()
end

