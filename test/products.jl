# Minimal downstream model implementing the reader contract. Its extra
# type needs no coupling to the storage adapter or product implementation.
struct ProductTestSelector
    source::Aegis.FieldRef
    centering::Symbol
    name::String
    units::String
end
function Aegis.resolve_field(s::Aegis.Snapshot, f::ProductTestSelector)
    resolve_field(s,f.source)
    f
end
function Aegis.read_field(s::Aegis.Snapshot, f::ProductTestSelector, b::Integer; ranges)
    2 .* read_field(s,f.source,b; ranges)
end

# Two native layers along x. Independent of the production index mapping.
function layered_product_fixture(path, order)
    nx = 2 ÷ length(order)
    write_fixture(path; order, blocksize=(nx,2,1), fields=("rho","q"))
    h5open(path,"r+") do h
        HDF5.delete_attribute(h,"RootGridX1")
        attributes(h)["RootGridX1"] = [0.,4.,-1.]
        for (b,loc) in pairs(order)
            r = (loc-1)*nx+1:loc*nx
            faces = [0.,1.,4.][first(r):last(r)+1]
            h["x1f"][:,b] = faces
            h["x1v"][:,b] = (faces[1:end-1]+faces[2:end])/2
            for (d,field) in (("gas",[3,1]),("extras",[2,6]))
                slab = repeat(reshape(Float32.(field[r]),:,1,1),1,2,1)
                h[d][:,:,:,b,1] = slab
            end
        end
    end
    path
end

@testset "native distributions and spatial products" begin
    mktempdir() do dir
        for order in ([1],[2,1])
            path = layered_product_fixture(joinpath(dir,"layers.h5"),order)
            open_snapshot(path) do s
                mesh = CartesianMesh(s)
                q = interpret_field(s,"q";units="u",alias="scalar")
                custom = ProductTestSelector(resolve_field(s,q), :cell, "twice-q", "2q unit")
                custom_plane = project(mesh,custom;axis=1)
                @test all(==(40.),custom_plane.values)
                @test custom_plane.metadata.field.units == "2q unit"
                @test custom_plane.metadata.field.kind === :user_defined
                x = BinAxis(q,[1.,2.,4.,6.];measure=:dex)
                y = BinAxis("rho",[0.,2.,3.])
                dist = distribution(mesh,x,y;weight="rho")
                # Cross section has area 3, each layer has equal column mass.
                @test dist.totals == [0 0; 0 9; 9 0]
                dx = log10.([2.,2.,1.5])
                dy = [2.,1.]
                @test dist.widths[1] ≈ dx
                @test dist.pdf ≈ dist.totals ./ (18 .* dx .* dy')
                @test sum(dist.pdf .* dx .* dy') ≈ 1
                @test dist.included_weight == 18
                @test dist.excluded_weight == 0
                @test dist.metadata.bin_fields[1].units == "u"
                @test distribution(mesh,x;weight=:count).totals == [0,2,2]
                narrow = distribution(mesh,BinAxis("q",[6.,7.]),BinAxis("rho",[1.,2.]))
                @test narrow.excluded_weight == 3  # outside both axes, counted once
                @test narrow.included_weight == 9
                a = BinAxis("q",[2.,6.,7.])
                ratio = cell_map(/, s, q, "rho";name="q/rho")
                binned = conditional_stats(mesh,a,ratio)
                @test binned.numerator == [2,54]
                @test binned.denominator == [3,9]
                ac = BinAxis("q",[1.,10.])
                @test only(conditional_stats(mesh,ac,ratio).mean) ≈ 14/3
                n = only(conditional_stats(mesh,ac,q).numerator)
                d = only(conditional_stats(mesh,ac,"rho").numerator)
                @test n/d ≈ 10/3
                @test summarize(mesh,ratio).integral == 56

                col = project(mesh,q;axis=1)
                @test col.faces == ([0,1,3],[0,1])
                @test col.values == fill(20.,2,1)
                @test col.path_length == 4
                @test project(mesh,q;axis=1,mode=:mean).values == fill(5.,2,1)
                m = project(mesh,q;axis=1,mode=:mean,weight="rho")
                @test m.denominator == fill(6.,2,1)
                @test m.values == fill(4.,2,1)
                @test sum(col.values .* diff(col.faces[1]) .* diff(col.faces[2])') == 60
                @test summarize(mesh,q).integral == 60
                @test project(mesh,q;axis=2).values == reshape([6.,18.],2,1)
                @test project(mesh,q;axis=3).values == [2 2;6 6]
                @test slice_plane(mesh,q;axis=1,at=1.).values == fill(6.,2,1)
                @test slice_plane(mesh,q;axis=1,at=prevfloat(1.)).values == fill(2.,2,1)
                @test slice_plane(mesh,q;axis=2,at=0.).values == reshape([2.,6.],2,1)
                @test slice_plane(mesh,q;axis=3,at=.3).values == [2 2;6 6]
                r = Region(x1=(1,4),x2=(1,3))
                @test slice_plane(mesh,q;axis=1,at=2.,region=r).values == fill(6.,1,1)
                @test project(mesh,q;axis=1,region=r).values == fill(18.,1,1)
                @test_throws ArgumentError slice_plane(mesh,q;axis=1,at=4.)
                @test_throws ArgumentError slice_plane(mesh,q;axis=1,at=.2,region=r)
                @test_throws ArgumentError project(mesh,q;axis=4)
                @test_throws ArgumentError project(mesh,q;axis=1,weight="rho")

                # A tiny bin probability underflows Float64 before width
                # normalization, yet its per-unit-x density is representable.
                tiny_coordinate = cell_map(t -> t==2 ? 0.0 : 1.0, s, q;name="test coordinate")
                tiny_weight = cell_map(t -> t==2 ? 1e-300 : 1e100, s,q;name="test weight")
                narrow_pdf = distribution(mesh,BinAxis(tiny_coordinate,[0.,1e-250,2.]);weight=tiny_weight)
                @test narrow_pdf.pdf[1] ≈ 1e-150/3

                emptybins = distribution(mesh,BinAxis("q",[7.,8.]))
                @test emptybins.excluded_weight == 12 && only(emptybins.totals)==0 && isnan(only(emptybins.pdf))
                emptycond = conditional_stats(mesh,BinAxis("q",[7.,8.]),"rho")
                @test isnan(only(emptycond.mean))
                maybezero = cell_map(x -> x-2,s,q;name="q-2")
                logpositive = distribution(mesh, BinAxis(maybezero,[2.,8.];measure=:dex))
                @test logpositive.excluded_weight == 3 && logpositive.included_weight == 9
                neg = cell_map(-,s,q;name="-q")
                @test_throws DomainError distribution(mesh,x;weight=neg)
                @test_throws DomainError distribution(mesh,BinAxis(cell_map(t->NaN,s,q;name="nan"),[0,1]))
                @test all(isnan,project(mesh,q;axis=2,mode=:mean,weight=cell_map(zero,s,q;name="zero")).values)

                # A derived selector uses the same cell-reading path and keeps
                # its EOS and role in the export.
                model = NewtonianMHD(; eos=Isothermal(2.))
                cs = quantity(s,model,:sound_speed)
                cd = distribution(mesh,BinAxis(cs,[1.,3.]))
                @test only(cd.totals) == 12
                @test cd.metadata.bin_fields[1].eos.sound_speed == 2
                for p in (dist,binned,m,slice_plane(mesh,ratio;axis=1,at=0.),cd)
                    mktempdir(dir) do output
                        file = export_product(joinpath(output,"data.h5"),p)
                        h5open(file,"r") do h
                            @test read(h["kind"]) == string(p.kind)
                            @test read(h["metadata/source_path"]) == s.path
                            @test read(h["metadata/region_bounds/1/1"]) == -Inf
                            @test HDF5.read_attribute(h,"aegis_product_schema")==1
                            key = hasproperty(p,:values) ? :values : hasproperty(p,:totals) ? :totals : :numerator
                            @test read(h[string(key)]) == getproperty(p,key)
                        end
                        @test_throws ArgumentError export_product(file,p)
                    end
                end
            end
        end
        @test_throws ArgumentError BinAxis("q",[0,0])
        @test_throws ArgumentError BinAxis("q",[0,1];measure=:dex)
        @test_throws ArgumentError BinAxis("q",[0,Inf])
    end
end
