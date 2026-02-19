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
    marked_flag[marked_elem] .= 1
    @show size(marked_elem)

    # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TETRA, cny[:, i]) for i = 1:size(cny, 2)]
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 1:size(cny, 2)]

    vtk_grid("result/vtu/mesh", coor, cells) do vtk
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

function write_mesh_deformation(data_dir)
    # cny = reshape(reinterpret(Int32, read(data_dir * "cny_quad.bin")), (11, :))
    # material = cny[11, :]
    # cny = cny[1:10, :]
    cny = reshape(reinterpret(Int32, read(data_dir * "cny_linear.bin")), (5, :))
    material = cny[5, :]
    cny = cny[1:4, :]
    # coor = reshape(reinterpret(Float64, read(data_dir * "coor_quad.bin")), (3, :))
    coor = reshape(reinterpret(Float64, read(data_dir * "coor_linear.bin")), (3, :))
    # eta = reshape(reinterpret(Float64, read(data_dir * "eta_quad.bin")), (1, :))
    # displacement = reshape(reinterpret(Float64, read(data_dir * "displacement_quad.bin")), (3, :))
    displacement = reshape(reinterpret(Float64, read(data_dir * "displacement_linear.bin")), (3, :))

    # marked_elem = reshape(reinterpret(Int32, read(data_dir * "marked_elem_quad.bin")), (1, :))
    # load_elem = reshape(reinterpret(Int32, read(data_dir * "load_elem.bin")), (1, :))
    surf_elem = reshape(reinterpret(Int32, read(data_dir * "surf_elem.bin")), (1, :))
    partition = reshape(reinterpret(Int32, read(data_dir * "partition.bin")), (1, :))
    @show size(coor)
    @show size(cny)
    # @show size(eta)
    @show size(displacement)

    @show sum(abs.(coor))

    # # -3 -> -(3^(0.5))
    # # -4 -> -(4^(0.5))
    # sign = displacement ./ abs.(displacement)
    # displacement = sign .* (abs.(displacement)).^(0.5)

    # coor .+= displacement .* 10^(7.0)

    @show sum(abs.(coor))

    @show "number of element: ", size(cny, 2)
    # marked_flag = zeros(Int32, size(cny, 2))
    # marked_flag[marked_elem] .= 1
    # @show size(marked_elem)

    # cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TETRA, cny[:, i]) for i = 1:size(cny, 2)]
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, cny[:, i]) for i = 1:size(cny, 2)]

    vtk_grid("result/vtu/mesh", coor, cells) do vtk
        # vtk["eta"] = eta[1, :]
        vtk["material"] = material
        vtk["displacement_x"] = displacement[1, :]
        vtk["displacement_y"] = displacement[2, :]
        vtk["displacement_z"] = displacement[3, :]
        # vtk["marked_flag"] = marked_flag
        # vtk["load_elem"] = load_elem
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
    original = reshape(reinterpret(Int32, read(data_dir * "new_original.bin")), (1, :))
    node_id = 0:size(coor, 2) - 1
    elem_id = 0:size(cny, 2) - 1
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


    vtk_grid("result/new_mesh", coor, cells) do vtk
        vtk["new_elem"] = new_elem
        vtk["material"] = material
        vtk["original"] = original
        vtk["node_id"] = node_id
        vtk["elem_id"] = elem_id
    end
end

function write_mesh_gflib(data_dir)
    for iserver = 1:4
        @show "server: ", iserver
        cny = reshape(reinterpret(Int32, read(data_dir * "cny_quad_0$(iserver).bin")), (11, :))
        material = cny[11, :]
        cny = cny[1:10, :]
        coor = reshape(reinterpret(Float64, read(data_dir * "coor_quad_0$(iserver).bin")), (3, :))
        @show "nnode: ", size(coor, 2)
        @show "nelem: ", size(cny, 2)

        cells = [MeshCell(VTKCellTypes.VTK_QUADRATIC_TETRA, cny[:, i]) for i = 1:size(cny, 2)]
        vtk_grid(data_dir * "mesh_0$(iserver)", coor, cells) do vtk
            vtk["material"] = material
        end
    end
end

write_mesh("result/")
# println("------------")
# write_new_mesh("result/")

# write_mesh_gflib("gflib/local/")

# write_mesh_deformation("result/")
