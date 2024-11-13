using Test
using AtomicSymmetries
using PyCall
include("define_cell.jl")

function test_load_spglib()
    positions, cell, types = get_pm3m_perovskite()

    # Get the symmetry group
    pm3m_group = get_symmetry_group_from_spglib(positions, cell, types)

    # Check that the symmetries are 48
    @test length(pm3m_group) == 48
end

function test_gold_unit_cell_spglib()
    pos, cell, types = get_gold_unit_cell()
    gold_group = get_symmetry_group_from_spglib(pos, cell, types)
    @test length(gold_group) == 12
end

function test_gold_supercell_spglib()
    pos, cell, types = get_gold_unit_cell()
    positions, cell, types = get_supercell(pos, cell, types, [2, 2, 2])

    gold_group = get_symmetry_group_from_spglib(positions, cell, types)
    println("Positions: ", positions')
    println("Cell: ", cell')
    println("Types: ", types)
    println("Gold group: ", length(gold_group))
    @test length(gold_group) == 12*8
end


function test_python_spglib()
    pos, cell, types = get_gold_unit_cell()
    gold_group = get_symmetry_group_from_spglib(pos, cell, types; spglib_py_module = pyimport("spglib"))

    @test length(gold_group) == 12
end


function test_symmetries_supercell()
    positions, cell, types = get_pm3m_supercell()

    @show positions'

    # Get the symmetry group
    pm3m_group_supercell = get_symmetry_group_from_spglib(positions, cell, types)

    @test length(pm3m_group_supercell) == 48*8

    nat = size(positions, 2)
    nat_uc = nat ÷ 8
    my_matrix = rand(3*nat, 3*nat)
    my_matrix .+= my_matrix'

    # Now symmetrize the dynamical matrix
    pm3m_group_supercell.symmetrize_fc!(my_matrix)

    # Print the IRT of the translations
    for i in 1:8
        index = 48*(i-1) + 1
        println("Symmetry: ")
        @show pm3m_group_supercell.symmetries[index]
        println("IRT of the translation $i: ", pm3m_group_supercell.irt[index])
    end

    # Check a simple translational symmetry
    @test abs(my_matrix[1,1] - my_matrix[1 + 3nat_uc, 1 + 3nat_uc]) < 1e-12
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_load_spglib()
    test_symmetries_supercell()
end
