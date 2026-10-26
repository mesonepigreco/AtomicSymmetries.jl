using AtomicSymmetries, LinearAlgebra, Random, Test
const AS = AtomicSymmetries
Random.seed!(992)
@testset "AtomicSymmetries neighboring checkout" begin
    @testset "Fourier API compatibility" begin
        q=zeros(3,1); R=zeros(3,1); itau=[1]
        @test applicable(AS.vector_r2q!,zeros(ComplexF64,3,1),zeros(3),q,itau,R)
        @test applicable(AS.vector_r2q!,zeros(ComplexF64,3,1),zeros(3),q,itau,R)
    end
    @testset "Hermitian projection" begin
        M = zeros(ComplexF64,3,3,1); M[:,:,1] .= (1+2im).*Matrix(I,3,3)
        AS.impose_hermitianity_q!(M,[1])
        println("Hermitian projection imaginary diagonal: ",imag(M[1,1,1]))
        @test M[:,:,1] ≈ M[:,:,1]'
    end
    @testset "Matrix symmetrization alias contract" begin
        group=AS.get_identity_symmetry_group(Float64;dims=3,n_atoms=1)
        qs=AS.SymmetriesQSpace(group,zeros(3,1))
        M=reshape(ComplexF64.(Matrix(I,3,3)),3,3,1)
        original=copy(M)
        AS.symmetrize_matrix_q!(M,M,qs)
        println("Identity symmetry with aliased input norm: ",norm(M))
        @test M ≈ original
    end
    @testset "Single operation vector and stress covariance" begin
        group=AS.get_identity_symmetry_group(Float64;dims=3,n_atoms=1)
        S=[0. -1 0; 1 0 0; 0 0 1]
        AS.add_symmetry!(group,S;irt=[1])
        cell=Matrix{Float64}(I,3,3)
        recip=zeros(3,3); AS.get_reciprocal_lattice!(recip,cell)
        v=[1.,2,3]; out=zeros(3)
        AS.rotate_vector!(out,v,cell,recip,group,2)
        stress=v*v'; rotated=zeros(3,3)
        AS.rotate_matrix!(rotated,stress,cell,recip,group,2)
        @test out ≈ S*v
        println("Stress covariance error: ",maximum(abs.(rotated-out*out')))
        @test rotated ≈ out*out'
    end
end
@testset "Coset decomposition with atom-resolved perturbation" begin
    positions=[0. 0.25 0.5 0.75; 0 0 0 0; 0 0 0 0]
    cell=Matrix(Diagonal([1.,1.3,1.7]))
    group=AS.get_symmetry_group_from_spglib(positions,cell,[1,1,1,1])
    result = try
        AS.CosetSymmetryGroup(group,[1.,0,0,2,0,0,3,0,0,4,0,0],cell)
    catch e
        e
    end
    println("Atom-resolved coset result: ",result isa Exception ? sprint(showerror,result) : length(result))
    @test !(result isa Exception)
end

@testset "Lattice gauge fractional translation covariance" begin
    positions=[0. 0.5; 0 0; 0 0]
    cell=Matrix{Float64}(I,3,3); recip=zeros(3,3)
    get_reciprocal_lattice!(recip,cell)
    group=get_symmetry_group_from_spglib(positions,cell,[1,1])
    k=findfirst(i -> group.symmetries[i] ≈ cell && group.irt[i]==[2,1],1:length(group))
    @test k !== nothing
    q=[0. 1/3 2/3;0 0 0;0 0 0]; R=[0. 0 1 1 2 2;zeros(2,6)]
    itau=repeat([1,2],3); qs=SymmetriesQSpace(group,q)
    # Half-cell translation: each source atom advances by one site in a six-site ring.
    permutation=[2,3,4,5,6,1]
    u=randn(18); rotated=zeros(18)
    for a in 1:6
        rotated[3permutation[a]-2:3permutation[a]] .= u[3a-2:3a]
    end
    uq=zeros(ComplexF64,6,3); expected=similar(uq); actual=similar(uq)
    vector_r2q!(uq,u,q,itau,R); vector_r2q!(expected,rotated,q,itau,R)
    rotate_vector!(actual,uq,cell,recip,qs,k)
    @test actual ≈ expected
    # Independent block-circulant matrix, and its real-space pullback P' M P.
    blocks=randn(6,6,3); M=zeros(18,18)
    for a in 0:2,b in 0:2
        M[6a+1:6a+6,6b+1:6b+6] .= blocks[:,:,mod(b-a,3)+1]
    end
    M .= (M+M')/2
    indices=vcat([collect(3a-2:3a) for a in permutation]...)
    Mq=zeros(ComplexF64,6,6,3); out=similar(Mq); target=similar(Mq)
    matrix_r2q!(Mq,M,q,itau,R)
    matrix_r2q!(target,M[indices,indices],q,itau,R)
    fill!(out,0)
    apply_symmetry_matrixq!(out,Mq,group.symmetries[k],group.irt[k],qs.irt_q[k],group.unit_cell_translations[k],q)
    @test out ≈ target
    rebuilt=deepcopy(group)
    AS.initialize_uc_translations!(rebuilt,positions)
    @test rebuilt.unit_cell_translations[k] ≈ group.unit_cell_translations[k]
end

@testset "Hermitian and time reversal independent projections" begin
    q=[0. 1/3 2/3;0 0 0;0 0 0]; minus=[1,3,2]
    M=randn(ComplexF64,6,6,3)
    impose_hermitianity_q!(M,minus)
    for iq in 1:3
        @test M[:,:,iq] ≈ M[:,:,iq]'
        @test M[:,:,iq] ≈ conj.(M[:,:,minus[iq]])
    end
    projected=copy(M); impose_hermitianity_q!(M,minus)
    @test M ≈ projected
end
