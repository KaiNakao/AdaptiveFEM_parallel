subroutine read_nodal_force(myid, nnode, nelem, ni, coor, cny, rv)
    implicit none
    integer, intent(in) :: myid, nnode, nelem, ni, cny(10, nelem)
    double precision, intent(in) :: coor(3, nnode)
    double precision, intent(inout) :: rv(3*(nnode+ni))

    integer :: nnode_tmp, inode, node_id
    double precision :: fx, fy, fz

    character(len=100) :: filename
    
    rv = 0d0
    write(filename, '("mdata/", I6.6, ".nodal_force.dat")') myid
    open(10, file=filename, status="old")
    read(10, *) ! number of nodes 
    read(10, *) nnode_tmp
    read(10, *) ! node_id fx fy fz
    if (nnode_tmp > 0) then
        print *, "read nodal force"
    end if
    do inode = 1, nnode_tmp
        read(10, *) node_id, fx, fy, fz
        print *, node_id, fx, fy, fz
        rv(3*node_id - 2) = fx
        rv(3*node_id - 1) = fy
        rv(3*node_id) = fz
    end do
    close(10)

end subroutine read_nodal_force