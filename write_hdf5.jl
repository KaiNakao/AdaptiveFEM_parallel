using HDF5

function main(data_dir)
    # cny = reshape(reinterpret(Int32, read("/home/nakao/seismic_wave/input/cny_linear.bin")), (4, :))
    # cny = Int64.(cny)

    # cny_with_material = zeros(Int64, 5, size(cny, 2))
    # cny_with_material[1:4, :] = copy(cny)
    # cny_with_material[5, :] = ones(Int64, size(cny, 2))
    cny_with_material = reshape(reinterpret(Int32, read(data_dir * "new_connectivity.bin")), (5, :))
    cny_with_material = Int64.(cny_with_material)

    coor = reshape(reinterpret(Float64, read(data_dir * "new_coordinates.bin")), (3, :))
    @show size(cny_with_material)
    @show size(coor)

    h5open(data_dir * "tet4vox8model.h5", "w") do file
        file["/coor"] = reshape(coor, 3*size(coor, 2))
        file["/conn_tet4"] = reshape(cny_with_material, 5*size(cny_with_material, 2))
    end

    open(data_dir * "tet4vox8setting.dat", "w") do file
        write(file, "num of node\n")
        write(file, string(size(coor, 2)))
        write(file, "\n")
        write(file, "num of tet elem\n")
        write(file, string(size(cny_with_material, 2)))
        write(file, "\n")
        write(file, "num of vox elem\n")
        write(file, string(0))
        write(file, "\n")
        write(file, "num of material properties\n")
        write(file, string(maximum(cny_with_material[5, :])))
    end
end 

main("/data6/itou/AFEM/data/analysis_result_iburi_large/")