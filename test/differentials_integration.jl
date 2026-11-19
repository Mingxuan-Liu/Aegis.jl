# A direct two-secant reconstruction at a real block seam, independent of the
# block gatherer and quadratic implementation. This file can run standalone.
using Aegis, Test
@testset "vorticity at a filament-run MeshBlock seam" begin
    root = joinpath(@__DIR__,"test_data","3D_Athena","128_eta5e-3")
    input = read_athinput(only(filter(f -> startswith(basename(f),"athinput"),readdir(root;join=true))))
    @test all(parameter(input,"mesh","$(side)x$(d)_bc") == "periodic" for side in ("i","o"), d in 1:3)
    record = last(only(values(discover_snapshots(root))))
    open_snapshot(record.path) do s
        mesh = CartesianMesh(s)
        grid = DifferentialGrid(mesh;boundary=(:periodic,:periodic,:periodic))
        omega = curl(grid,("vel1","vel2","vel3"))
        b = only(b for b in 1:s.nblocks if s.locations[:,b]==[1,0,0])
        left = only(b for b in 1:s.nblocks if s.locations[:,b]==[0,0,0])
        j,k = 12,17
        sample(field,block,i,j,k) = Float64(only(read_field(s,field,block;ranges=(i:i,j:j,k:k))))
        cx,cy,_ = block_coordinates(s,b).centers
        lx = block_coordinates(s,left).centers[1]
        dy_dx = (sample("vel2",b,2,j,k)-sample("vel2",left,s.blocksize[1],j,k))/(cx[2]-lx[end])
        dx_dy = (sample("vel1",b,1,j+1,k)-sample("vel1",b,1,j-1,k))/(cy[j+1]-cy[j-1])
        measured = only(read_field(s,omega[3],b;ranges=(1:1,j:j,k:k)))
        @test measured ≈ dy_dx-dx_dy atol=8eps(Float64)*max(abs(dy_dx),abs(dx_dy))
        println("t=",s.time,", logical block (1,0,0), native x1 seam: omega3=", measured)
    end
end
