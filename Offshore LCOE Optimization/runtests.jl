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
    @test calculate_50year_wind_load(9.74, 125, 11.4) == 2087529.8529107259
    @test calculate_50year_wind_load(7, 125, 11.4) == 1618939.5140525973
    @test calculate_50year_wind_load(7, 100, 11.4) == 1070554.8426587582
    @test calculate_50year_wind_load(7, 125, 10) == 1590230.7187345193
end

@testset "50 year Wind Moment" begin
    @test calculate_50year_wind_moment(9.74, 23.5, 125, 90, 11.4) == 799654404.2806149
    @test calculate_50year_wind_moment(7, 23.5, 125, 90, 11.4) == 620155017.6017731
    @test calculate_50year_wind_moment(8, 26, 125, 90, 10) == 691345156.9397321
    @test calculate_50year_wind_moment(8, 30, 125, 80, 11.4) == 661962945.240009
end

@testset "Monopile Diameter" begin
    @test isapprox(monopile_diameter(355000000, 1.1, 799654404.2806149), 6.677656219964458; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 620155017.6017731), 6.119550254228244; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 691345156.9397321), 6.352343360869822; atol = 0.001)
    @test isapprox(monopile_diameter(355000000, 1.1, 661962945.240009), 6.25829318500689; atol = 0.001)
end

@testset "Monopile Thickness" begin
    @test isapprox(monopile_thickness(6.677656219964458), 0.07312656219964457; 0.001)
    @test isapprox(monopile_thickness(5), 0.056350000000000004; 0.001)
    @test isapprox(monopile_thickness(7), 0.07635; 0.001)
    @test isapprox(monopile_thickness(6.119550254228244), 0.06754550254228243; 0.001)
end

@testset "Monopile Moment" begin
    @test isapprox(calculate_pile_moment(6.677656219964458, 0.07312656219964457), 8.27295623552583; 0.001)
    @test isapprox(calculate_pile_moment(5, 0.056350000000000004), 2.6736032113211996; 0.001)
    @test isapprox(calculate_pile_moment(7, 0.07635), 9.95117225211555; 0.001)
    @test isapprox(calculate_pile_moment(6.119550254228244, 0.06754550254228243), 5.879685598856954; 0.001)
end

@testset "Embedment Length" begin
    @test isapprox(pile_embedment_length(8.27295623552583), 26.56783356260234; 0.001)
    @test isapprox(pile_embedment_length(2.6736032113211996), 21.195486401869456; 0.001)
    @test isapprox(pile_embedment_length(9.95117225211555), 27.567592805981995; 0.001)
    @test isapprox(pile_embedment_length(5.879685598856954), 24.813887979958285; 0.001)
end

@testset "Monopile Mass" begin
    @test isapprox(monopile_mass(6.677656219964458, 0.07312656219964457, (23.5+10+26.56783356260234)), 397.01164557176816; 0.001)
    @test isapprox(monopile_mass(6.119550254228244, 0.06754550254228243, (23.5+10+24.813887979958285)), 326.2353875559374; 0.001)
end