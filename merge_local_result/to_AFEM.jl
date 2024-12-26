using HDF5
using LinearAlgebra

work_dir = "./"
nproc = 16

neighbor_arr = []
nnode_arr = []
nelem_arr = []
cny_quad_arr = []
coor_quad_arr = []
displacement_quad_arr = []
load_elem_arr = []
surf_elem_arr = []

for iproc in 0:nproc - 1
    filename = work_dir * "mdata/" * lpad(iproc, 6, "0") * ".data.h5"
    fid = h5open(filename, "r")

    dat = read(fid["MPInode"])
    n = dat[1]
    rank = zeros(Int32, n)
    len = zeros(Int32, n)
    neighbor = Dict{Int32, Vector{Int32}}()
    for i in 1:n
        rank[i] = dat[2*i]
        len[i] = dat[2*i + 1]
    end
    pt = 2*n + 2
    for i in 1:n
        neighbor[rank[i]] = dat[pt:pt + len[i] - 1]
        pt += len[i]
    end
    push!(neighbor_arr, neighbor)

    dat = read(fid["setting"])
    nnode = dat[1]
    nelem = dat[2]
    push!(nnode_arr, nnode)
    push!(nelem_arr, nelem)

    dat = read(fid["conn"])
    cny_quad = reshape(dat, (11, :))
    push!(cny_quad_arr, cny_quad)

    dat = read(fid["coor"])
    coor_quad = reshape(dat, (3, :))
    push!(coor_quad_arr, coor_quad)

    close(fid)

    filename = work_dir * "displacement/" * lpad(iproc, 6, "0") * ".bin"
    displacement_quad = reshape(reinterpret(Float64, read(filename)), (3, :))
    push!(displacement_quad_arr, displacement_quad)

    filename = work_dir * "mdata/" * lpad(iproc, 6, "0") * ".load_elem.bin"
    load_elem = reinterpret(Int32, read(filename))
    push!(load_elem_arr, load_elem)
    for ielem = 1:nelem
        if (load_elem[ielem] == 1)
            @show "load_elem", iproc, ielem
        end
    end

    filename = work_dir * "surf_mesh/" * lpad(iproc, 6, "0") * "_cny.bin"
    tmp = reshape(reinterpret(Int32, read(filename)), (4, :))
    surf_elem = zeros(Int32, nelem)
    surf_elem[tmp[4, :]] .= 1
    push!(surf_elem_arr, surf_elem)
end

rank_added = []
local_to_global_arr = []
node_id_cnt = 0
for iproc in 0:nproc - 1
    @show iproc, rank_added, keys(neighbor_arr[iproc + 1])

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
            # find inode in neighbor[jproc]
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
            global node_id_cnt
            node_id_cnt += 1
            local_to_global[inode] = node_id_cnt
        end
    end
    push!(rank_added, iproc)
    push!(local_to_global_arr, local_to_global)
end

# reorder nodes
vertex_nodes = Set()
edge_nodes = Set()
for iproc in 0:nproc - 1
    local_to_global = local_to_global_arr[iproc + 1]
    cny_quad = cny_quad_arr[iproc + 1]
    for ielem in 1:nelem_arr[iproc + 1]
        for i in 1:4
            push!(vertex_nodes, local_to_global[cny_quad[i, ielem]])
        end
        for i in 5:10
            push!(edge_nodes, local_to_global[cny_quad[i, ielem]])
        end
    end
end
vertex_nodes_arr = sort(collect(vertex_nodes))
edge_nodes_arr = sort(collect(edge_nodes))

old_to_new = Dict{Int32, Int32}()
for (i, inode) in enumerate(vertex_nodes_arr)
    old_to_new[inode] = i
end
for (i, inode) in enumerate(edge_nodes_arr)
    old_to_new[inode] = i + length(vertex_nodes_arr)
end

for iproc in 0:nproc - 1
    local_to_global = local_to_global_arr[iproc + 1]
    tmp = Dict{Int32, Int32}()
    for inode in keys(local_to_global)
        tmp[inode] = old_to_new[local_to_global[inode]]
    end
    local_to_global_arr[iproc + 1] = tmp
end

cny_quad_global = zeros(Int32, (11, sum(nelem_arr)))
coor_quad_global = zeros(Float64, (3, node_id_cnt))
displacement_quad_global = zeros(Float64, (3, node_id_cnt))
load_elem_global = zeros(Int32, sum(nelem_arr))
surf_elem_global = zeros(Int32, sum(nelem_arr))
partition = zeros(Int32, sum(nelem_arr))
pt = 0
for iproc in 0:nproc - 1
    global pt
    nelem = nelem_arr[iproc + 1]
    nnode = nnode_arr[iproc + 1]
    local_to_global = local_to_global_arr[iproc + 1]
    cny_quad = cny_quad_arr[iproc + 1]
    coor_quad = coor_quad_arr[iproc + 1]
    load_elem = load_elem_arr[iproc + 1]
    surf_elem = surf_elem_arr[iproc + 1]
    displacement_quad = displacement_quad_arr[iproc + 1]

    for ielem in 1:nelem
        for inode in 1:10
            cny_quad_global[inode, pt + ielem] = local_to_global[cny_quad[inode, ielem]]
        end
        cny_quad_global[11, pt + ielem] = cny_quad[11, ielem]
        partition[pt + ielem] = iproc
        load_elem_global[pt + ielem] = load_elem[ielem]
        surf_elem_global[pt + ielem] = surf_elem[ielem]
        if (pt + ielem == 143713)
            @show "surf_elem", iproc, ielem
        end
    end
    pt += nelem

    for inode in 1:nnode
        coor_quad_global[:, local_to_global[inode]] = coor_quad[:, inode]
        displacement_quad_global[:, local_to_global[inode]] = displacement_quad[:, inode]
    end 
end 

coor_linear_global = zeros(Float64, (3, length(vertex_nodes_arr)))
for inode in 1:length(vertex_nodes_arr)
    coor_linear_global[:, inode] = coor_quad_global[:, inode]
end
cny_linear_global = zeros(Int32, (5, sum(nelem_arr)))
for ielem in 1:sum(nelem_arr)
    for i in 1:4
        cny_linear_global[i, ielem] = cny_quad_global[i, ielem]
    end
    cny_linear_global[5, ielem] = cny_quad_global[11, ielem]
end

result_dir = work_dir * "result/"
if !isdir(result_dir)
    mkpath(result_dir)
end

open(work_dir * "result/coor_quad.bin", "w") do io
    write(io, coor_quad_global)
end
open(work_dir * "result/cny_quad.bin", "w") do io
    write(io, cny_quad_global)
end
open(work_dir * "result/partition.bin", "w") do io
    write(io, partition)
end
open(work_dir * "result/displacement_quad.bin", "w") do io
    write(io, displacement_quad_global)
end
open(work_dir * "result/coor_linear.bin", "w") do io
    write(io, coor_linear_global)
end
open(work_dir * "result/cny_linear.bin", "w") do io
    write(io, cny_linear_global)
end
open(work_dir * "result/load_elem.bin", "w") do io
    write(io, load_elem_global)
end
open(work_dir * "result/surf_elem.bin", "w") do io
    write(io, surf_elem_global)
end

nmaterial = maximum(cny_linear_global[5, :]) 
material = zeros(Float64, (3, nmaterial))
io = open(work_dir * "data/material.dat", "r")
lines = readlines(io)
for imaterial in 1:nmaterial
    material[1, imaterial] = parse(Float64, lines[imaterial * 11 - 9])
    material[2, imaterial] = parse(Float64, lines[imaterial * 11 - 8])
    material[3, imaterial] = parse(Float64, lines[imaterial * 11 - 7])
end
close(io)

open(work_dir * "result/material.bin", "w") do io
    write(io, material)
end

open(work_dir * "result/shape.dat", "w") do io
    write(io, "nelem\n")
    write(io, string(size(cny_linear_global, 2)))
    write(io, "\n")

    write(io, "nnode_linear\n")
    write(io, string(size(coor_linear_global, 2)))
    write(io, "\n")

    write(io, "nnode_quad\n")
    write(io, string(size(coor_quad_global, 2)))
    write(io, "\n")
    
    write(io, "nmaterial\n")    
    write(io, string(nmaterial))
    write(io, "\n")
end
