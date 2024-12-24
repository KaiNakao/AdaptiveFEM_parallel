using WriteVTK
function write_mesh(data_dir)
    # cny = reshape(reinterpret(Int32, read(data_dir * "cny_quad.bin")), (11, :))
    # material = cny[11, :]
    # cny = cny[1:10, :]
    cny = reshape(reinterpret(Int32, read(data_dir * "cny_linear.bin")), (5, :))
    material = cny[5, :]
    cny = cny[1:4, :]
    # coor = reshape(reinterpret(Float64, read(data_dir * "coor_quad.bin")), (3, :))
    coor = reshape(reinterpret(Float64, read(data_dir * "coor_linear.bin")), (3, :))
    eta = reshape(reinterpret(Float64, read(data_dir * "eta_quad.bin")), (1, :))
    # displacement = reshape(reinterpret(Float64, read(data_dir * "displacement_quad.bin")), (3, :))
    marked_elem = reshape(reinterpret(Int32, read(data_dir * "marked_elem_quad.bin")), (1, :))
    load_elem = reshape(reinterpret(Int32, read(data_dir * "load_elem.bin")), (1, :))
    surf_elem = reshape(reinterpret(Int32, read(data_dir * "surf_elem.bin")), (1, :))
    partition = reshape(reinterpret(Int32, read(data_dir * "partition.bin")), (1, :))
    @show size(coor)
    @show size(cny)
    # @show size(eta)
    # @show size(displacement)

    @show "number of element: ", size(cny, 2)

    marked_flag = zeros(Int32, size(cny, 2))
    marked_flag .= 0
    marked_flag[marked_elem] .= 1


    # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TETRA, cny[:, i]) for i = 1:size(cny, 2)]
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 1:size(cny, 2)]

    vtk_grid("./mesh", coor, cells) do vtk
        vtk["eta"] = eta[1, :]
        vtk["material"] = material
        # vtk["displacement_x"] = displacement[1, :]
        # vtk["displacement_y"] = displacement[2, :]
        # vtk["displacement_z"] = displacement[3, :]
        vtk["marked_flag"] = marked_flag
        vtk["load_elem"] = load_elem
        vtk["surf_elem"] = surf_elem
        vtk["partition"] = partition
        vtk["elem_id"] = 1:size(cny, 2)
    end
end 

function write_new_mesh(data_dir)
    cny = reshape(reinterpret(Int32, read(data_dir * "new_connectivity.bin")), (5, :))
    material = cny[5, :]
    cny = cny[1:4, :]
    coor = reshape(reinterpret(Float64, read(data_dir * "new_coordinates.bin")), (3, :))
    # eta = reshape(reinterpret(Float64, read(data_dir * "eta_quad.bin")), (1, :))
    # displacement = reshape(reinterpret(Float64, read(data_dir * "displacement_quad.bin")), (3, :))
    # marked_elem = reshape(reinterpret(Int32, read(data_dir * "marked_elem_quad.bin")), (1, :))
    # load_elem = reshape(reinterpret(Int32, read(data_dir * "load_elem.bin")), (1, :))
    # surf_elem = reshape(reinterpret(Int32, read(data_dir * "surf_elem.bin")), (1, :))
    # partition = reshape(reinterpret(Int32, read(data_dir * "partition.bin")), (1, :))
    # @show size(coor)
    # @show size(cny)
    # @show size(eta)
    # @show size(displacement)
    @show "number of element: ", size(cny, 2)
    @show "number of nodes: ", size(coor, 2)
    
    new_elem = zeros(Int32, size(cny, 2))
    for ielem = 1:size(cny, 2)
        for inode = 1:4
            if cny[inode, ielem] > 241032
                new_elem[ielem] = 1
                break
            end
        end
    end

    # marked_flag = zeros(Int32, size(cny, 2))
    # marked_flag .= 0
    # marked_flag[marked_elem] .= 1

    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 1:size(cny, 2)]


    vtk_grid("./new_mesh", coor, cells) do vtk
        vtk["new_elem"] = new_elem
        vtk["material"] = material
    end
end

function write_tmp()
    cny = zeros(Int32, (4, 3))
    cny[:, 1] .= [13196, 13275, 13274, 241092] .+ 1
    cny[:, 2] .= [13275, 13285, 13274, 241092] .+ 1
    cny[:, 3] .= [13285, 13196, 13274, 241092] .+ 1
    coor = zeros(Float64, (3, 241093))
    coor[:, 1 + 13196] .= [140000, 180000, 223505]
    coor[:, 1 + 13275] .= [130000, 190000, 215016]
    coor[:, 1 + 13274] .= [130000, 190000, 220000]
    coor[:, 1 + 13285] .= [140000, 190000, 222500]
    coor[:, 1 + 241092] .= [138845, 178845, 220816]

    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 1:1]
    vtk_grid("./tmp1", coor, cells) do vtk
    end
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 2:2]
    vtk_grid("./tmp2", coor, cells) do vtk
    end
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 3:3]
    vtk_grid("./tmp3", coor, cells) do vtk
    end
end

write_mesh("/data6/itou/AFEM/data/analysis_result_iburi_large/")
println("------------")
write_new_mesh("/data6/itou/AFEM/data/analysis_result_iburi_large/")
# write_tmp()