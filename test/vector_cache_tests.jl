struct TypeWithTrivialConstructor
    Int64::Int64
    TypeWithTrivialConstructor() = new()
    TypeWithTrivialConstructor(v) = new(v)
end

for type_k = [Int64, TypeWithTrivialConstructor]
    @testset "$(string(type_k))" begin
        vc = vectorCache{type_k}()
        @test length(vc) == -9999
        empty!(vc)
        @test isempty(vc)
        @test length(vc) == 0
        @test returnNext(vc) isa type_k

        empty!(vc)
        ind_max_1 = vc.ind_max
        for k = 1:ind_max_1
            if type_k == Int64
                addCacheItem!(vc, k)
                @test vc[k] == k
            elseif type_k == TypeWithTrivialConstructor
                addCacheItem!(vc, TypeWithTrivialConstructor(k))
                @test vc[k] == TypeWithTrivialConstructor(k)
            end
        end
        @test !isempty(vc)

        expand!(vc)
        ind_max_2 = vc.ind_max
        @test length(vc) == vc.ind_fill
        @test ind_max_1 < ind_max_2

        vc.ind_fill = vc.ind_max
        i_next = returnNext(vc)
        ind_max_3 = vc.ind_max
        @test ind_max_2 < ind_max_3
    end
end

@testset "No Trivial Constructor" begin
    @test_throws Exception vectorCache{Any}()
end
