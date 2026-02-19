using HDF5
using WriteVTK

function read_neighbor(dat::Vector{Int32})
    n = dat[1]
    ranks = Vector{Int32}(undef, n)
    lens = Vector{Int32}(undef, n)
    neighbor = Dict{Int32, Vector{Int32}}()
    for i in 1:n
        ranks[i] = dat[2 * i]
        lens[i] = dat[2 * i + 1]
    end
    pt = 2 * n + 2
    for i in 1:n
        neighbor[ranks[i]] = dat[pt:pt + lens[i] - 1]
        pt += lens[i]
    end
    return neighbor
end

function merge_rdata(work_dir::String, nproc::Int)
    neighbor_arr = Vector{Dict{Int32, Vector{Int32}}}()
    nnode_arr = Int[]
    nelem_arr = Int[]
    conn_arr = Vector{Matrix{Int32}}()
    coor_arr = Vector{Matrix{Float64}}()

    for iproc in 0:nproc - 1
        filename = joinpath(work_dir, "rdata", lpad(iproc, 6, "0") * ".data.h5")
        fid = h5open(filename, "r")

        neighbor = read_neighbor(read(fid["MPInode"]))
        push!(neighbor_arr, neighbor)

        conn = reshape(read(fid["conn"]), (5, :))
        coor = reshape(read(fid["coor"]), (3, :))

        nnode = size(coor, 2)
        nelem = size(conn, 2)
        push!(nnode_arr, nnode)
        push!(nelem_arr, nelem)
        push!(conn_arr, conn)
        push!(coor_arr, coor)
        close(fid)
    end

    rank_added = Int[]
    local_to_global_arr = Vector{Dict{Int32, Int32}}()
    node_id_cnt = 0
    for iproc in 0:nproc - 1
        local_to_global = Dict{Int32, Int32}()
        neighbor = neighbor_arr[iproc + 1]
        nnode = nnode_arr[iproc + 1]
        for inode in 1:nnode
            boundary = false
            b_rank = -1
            loc = nothing
            for jproc in keys(neighbor)
                if !(jproc in rank_added)
                    continue
                end
                loc = findfirst(x -> x == inode, neighbor[jproc])
                if loc != nothing
                    boundary = true
                    b_rank = jproc
                    break
                end
            end

            if boundary
                b_inode = neighbor_arr[b_rank + 1][iproc][loc]
                local_to_global[inode] = local_to_global_arr[b_rank + 1][b_inode]
            else
                node_id_cnt += 1
                local_to_global[inode] = node_id_cnt
            end
        end
        push!(rank_added, iproc)
        push!(local_to_global_arr, local_to_global)
    end

    # Compute max node id from actual mapping values
    max_id = 0
    for mp in local_to_global_arr
        for v in values(mp)
            if v > max_id
                max_id = v
            end
        end
    end
    node_id_cnt = max_id

    total_elems = sum(nelem_arr)
    conn_global = zeros(Int32, (5, total_elems))
    coor_global = zeros(Float64, (3, node_id_cnt))
    partition = zeros(Int32, total_elems)

    pt = 0
    for iproc in 0:nproc - 1
        nelem = nelem_arr[iproc + 1]
        nnode = nnode_arr[iproc + 1]
        local_to_global = local_to_global_arr[iproc + 1]
        conn = conn_arr[iproc + 1]
        coor = coor_arr[iproc + 1]

        for ielem in 1:nelem
            for inode in 1:4
                conn_global[inode, pt + ielem] = local_to_global[conn[inode, ielem]]
            end
            conn_global[5, pt + ielem] = conn[5, ielem]
            partition[pt + ielem] = iproc
        end
        pt += nelem

        for inode in 1:nnode
            idx = local_to_global[inode]
            if idx > size(coor_global, 2)
                new_global = zeros(Float64, (3, idx))
                new_global[:, 1:size(coor_global, 2)] = coor_global
                coor_global = new_global
            end
            coor_global[:, idx] = coor[:, inode]
        end
    end

    result_dir = joinpath(work_dir, "result")
    if !isdir(result_dir)
        mkpath(result_dir)
    end

    open(joinpath(result_dir, "coor_linear.bin"), "w") do io
        write(io, coor_global)
    end
    open(joinpath(result_dir, "cny_linear.bin"), "w") do io
        write(io, conn_global)
    end
    open(joinpath(result_dir, "partition.bin"), "w") do io
        write(io, partition)
    end

    vtu_dir = joinpath(result_dir, "vtu")
    if !isdir(vtu_dir)
        mkpath(vtu_dir)
    end

    cells = [MeshCell(VTKCellTypes.VTK_TETRA, conn_global[1:4, i]) for i in 1:total_elems]
    vtk_grid(joinpath(vtu_dir, "mesh"), coor_global, cells) do vtk
        vtk["material"] = conn_global[5, :]
        vtk["partition"] = partition
        vtk["elem_id"] = 1:total_elems
    end
end

function main()
    if length(ARGS) == 1
        work_dir = "./"
        nproc = parse(Int, ARGS[1])
        merge_rdata(work_dir, nproc)
    elseif length(ARGS) >= 2
        work_dir = ARGS[1]
        nproc = parse(Int, ARGS[2])
        merge_rdata(work_dir, nproc)
    else
        println("usage: julia merge_rdata_vtu.jl <nproc>")
        println("   or: julia merge_rdata_vtu.jl <work_dir> <nproc>")
    end
end

main()
