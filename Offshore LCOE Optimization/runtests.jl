#-------- UNIT TESTS for Monopile.jl--------#
#@test design_monopile(9.74, 23.5, 125, 90, 11.4)/3000 == 324
@testset "50 year extreme Wind Speed" begin
    @test calculate_50year_extreme_ws(9.74) == 37.43548973876968
    @test calculate_50year_extreme_ws(7) == 26.904356075091144
    @test calculate_50year_extreme_ws(8) == 30.74783551438988
end

@testset "50 year extreme Wind Gust" begin
    @test calculate_50year_extreme_gust(9.74, 125, 11.4) == 8.401646451820062
    @test calculate_50year_extreme_gust(7, 125, 11.4) == 6.038144267221809
    @test calculate_50year_extreme_gust(7, 100, 11.4) == 6.325538092410854
    @test calculate_50year_extreme_gust(7, 125, 10) == 6.038144267221809
end

@testset "Thrust Coefficient" begin
    @test calculate_thrust_coefficient(11.4) == 0.7082948599569098
    @test calculate_thrust_coefficient(7) == 1
    @test calculate_thrust_coefficient(8) == 1
    @test calculate_thrust_coefficient(10) == 0.8225
end

@testset "50 year Wind Load" begin
    @test calculate_50year_wind_load(9.74, 125, 11.4) == 72778861.9920054
    @test calculate_50year_wind_load(7, 125, 11.4) == 87627013.91109432
    @test calculate_50year_wind_load(7, 100, 11.4) == 30455906.948561136
    @test calculate_50year_wind_load(7, 125, 10) == 103988893.7694395
end

@testset "50 year Wind Moment" begin
    @test calculate_50year_wind_moment(9.74, 23.5, 125, 90, 11.4) == 799654404.2806149
    @test calculate_50year_wind_moment(7, 23.5, 125, 90, 11.4) == 620155017.6017731
    @test calculate_50year_wind_moment(8, 26, 125, 90, 10) == 691345156.9397321
    @test calculate_50year_wind_moment(8, 30, 125, 80, 11.4) == 661962945.240009
end