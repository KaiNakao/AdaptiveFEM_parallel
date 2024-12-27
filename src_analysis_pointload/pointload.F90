subroutine pointload(myid, nnode, coor, nelem, cny3d, ni, rv)
    use mpi
    implicit none
    
    integer, intent(in) :: myid, nnode, nelem, cny3d(10, nelem), ni
    double precision, intent(in) :: coor(3, nnode) 
    double precision, intent(inout) :: rv(3*(nnode + ni))
    integer :: ierr, surf_nelem, ielem, elem_id, inode, &
               found, found_cnt, load_elem, node_id(10), load_elem_arr(nelem)
    integer, allocatable :: surf_cny(:,:)
    double precision :: x, y, ex, ey, ez, xnode(3, 10), &
                        load_point(2), load_direction(3), &
                        dxdr(3, 3), drdx(3, 3), dx(3), detj, r1, r2, r3, &
                        dx3d(3, 2), load_point_proj(3), &
                        l0, l1, l2, l3, nvec(10), ftmp(30)
    character(len=100) :: filename

    if (myid == 0) then
        open(10, file="data/pointload.dat", status="old")
        read(10, *) 
        read(10, *) x, y, ex, ey, ez
        close(10)
    end if

    call mpi_bcast(x, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(y, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(ex, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(ey, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(ez, 1, mpi_double_precision, 0, mpi_comm_world, ierr)

    load_point(1) = x
    load_point(2) = y
    load_direction(1) = ex
    load_direction(2) = ey
    load_direction(3) = ez

    write(filename, '("surf_mesh/", I6.6, "_nelem.dat")') myid
    open(10, file=filename, status="old")
    read(10, *) surf_nelem
    close(10)

    allocate(surf_cny(4, surf_nelem))

    found = 0
    load_elem_arr = 0
    if (surf_nelem > 0) then
        write(filename, '("surf_mesh/", I6.6, "_cny.bin")') myid
        open(10, file=filename, status="old", access="stream", form="unformatted")
        read(10) surf_cny
        close(10)

        ! search elements that contains load point
        do ielem = 1, surf_nelem
            elem_id = surf_cny(4, ielem)
            node_id(1:3) = surf_cny(1:3, ielem)
            do inode = 1, 3
                xnode(:, inode) = coor(:, node_id(inode))
            end do

            dxdr(1, 1) = xnode(1, 2) - xnode(1, 1)
            dxdr(2, 1) = xnode(2, 2) - xnode(2, 1)
            dxdr(1, 2) = xnode(1, 3) - xnode(1, 1)
            dxdr(2, 2) = xnode(2, 3) - xnode(2, 1)

            detj = dxdr(1, 1) * dxdr(2, 2) - dxdr(1, 2) * dxdr(2, 1)

            drdx(1, 1) = dxdr(2, 2) / detj
            drdx(2, 1) = -dxdr(2, 1) / detj
            drdx(1, 2) = -dxdr(1, 2) / detj
            drdx(2, 2) = dxdr(1, 1) / detj

            dx(1) = load_point(1) - xnode(1, 1)
            dx(2) = load_point(2) - xnode(2, 1)

            r1 = drdx(1, 1) * dx(1) + drdx(1, 2) * dx(2)
            r2 = drdx(2, 1) * dx(1) + drdx(2, 2) * dx(2)

            if (-1d-8 <= r1 .and. r1 <= 1d0 + 1d-8 .and. -1d-8 <= r2 .and. r2 <= 1d0 + 1d-8.and. r1 + r2 <= 1d0 + 1d-8) then
                print *, "element found: ", elem_id
                if (found == 0) then
                    found = 1
                end if
                dx3d(1, 1) = xnode(1, 2) - xnode(1, 1)
                dx3d(2, 1) = xnode(2, 2) - xnode(2, 1)
                dx3d(3, 1) = xnode(3, 2) - xnode(3, 1)
                dx3d(1, 2) = xnode(1, 3) - xnode(1, 1)
                dx3d(2, 2) = xnode(2, 3) - xnode(2, 1)
                dx3d(3, 2) = xnode(3, 3) - xnode(3, 1)

                load_point_proj(1) = xnode(1, 1) + r1 * dx3d(1, 1) + r2 * dx3d(1, 2)
                load_point_proj(2) = xnode(2, 1) + r1 * dx3d(2, 1) + r2 * dx3d(2, 2)
                load_point_proj(3) = xnode(3, 1) + r1 * dx3d(3, 1) + r2 * dx3d(3, 2)

                print *, "load point projection: ", load_point_proj
                load_elem = elem_id

                load_elem_arr(elem_id) = 1
            end if
        end do
    end if

    rv = 0d0
    if (found == 0) then
    else
        node_id = cny3d(:, load_elem)
        do inode = 1, 10
            xnode(:, inode) = coor(:, node_id(inode))
        end do

        dxdr(1, 1) = xnode(1, 2) - xnode(1, 1)
        dxdr(2, 1) = xnode(2, 2) - xnode(2, 1)
        dxdr(3, 1) = xnode(3, 2) - xnode(3, 1)
        dxdr(1, 2) = xnode(1, 3) - xnode(1, 1)
        dxdr(2, 2) = xnode(2, 3) - xnode(2, 1)
        dxdr(3, 2) = xnode(3, 3) - xnode(3, 1)
        dxdr(1, 3) = xnode(1, 4) - xnode(1, 1)
        dxdr(2, 3) = xnode(2, 4) - xnode(2, 1)
        dxdr(3, 3) = xnode(3, 4) - xnode(3, 1)

        detj = dxdr(1, 1) * (dxdr(2, 2) * dxdr(3, 3) - dxdr(2, 3) * dxdr(3, 2)) &
             - dxdr(1, 2) * (dxdr(2, 1) * dxdr(3, 3) - dxdr(2, 3) * dxdr(3, 1)) &
             + dxdr(1, 3) * (dxdr(2, 1) * dxdr(3, 2) - dxdr(2, 2) * dxdr(3, 1))
        
        drdx(1, 1) = (dxdr(2, 2)*dxdr(3, 3) - dxdr(2, 3)*dxdr(3, 2)) / detj
        drdx(1, 2) = (dxdr(1, 3)*dxdr(3, 2) - dxdr(1, 2)*dxdr(3, 3)) / detj
        drdx(1, 3) = (dxdr(1, 2)*dxdr(2, 3) - dxdr(1, 3)*dxdr(2, 2)) / detj
        drdx(2, 1) = (dxdr(2, 3)*dxdr(3, 1) - dxdr(2, 1)*dxdr(3, 3)) / detj
        drdx(2, 2) = (dxdr(1, 1)*dxdr(3, 3) - dxdr(1, 3)*dxdr(3, 1)) / detj
        drdx(2, 3) = (dxdr(1, 3)*dxdr(2, 1) - dxdr(1, 1)*dxdr(2, 3)) / detj
        drdx(3, 1) = (dxdr(2, 1)*dxdr(3, 2) - dxdr(2, 2)*dxdr(3, 1)) / detj
        drdx(3, 2) = (dxdr(1, 2)*dxdr(3, 1) - dxdr(1, 1)*dxdr(3, 2)) / detj
        drdx(3, 3) = (dxdr(1, 1)*dxdr(2, 2) - dxdr(1, 2)*dxdr(2, 1)) / detj

        dx(1) = load_point_proj(1) - xnode(1, 1)
        dx(2) = load_point_proj(2) - xnode(2, 1)
        dx(3) = load_point_proj(3) - xnode(3, 1)

        r1 = drdx(1, 1) * dx(1) + drdx(1, 2) * dx(2) + drdx(1, 3) * dx(3)
        r2 = drdx(2, 1) * dx(1) + drdx(2, 2) * dx(2) + drdx(2, 3) * dx(3)
        r3 = drdx(3, 1) * dx(1) + drdx(3, 2) * dx(2) + drdx(3, 3) * dx(3)

        l0 = 1d0 - r1 - r2 - r3
        l1 = r1
        l2 = r2
        l3 = r3

        nvec(1) = l0*(2d0*l0 - 1d0);
        nvec(2) = l1*(2d0*l1 - 1d0);
        nvec(3) = l2*(2d0*l2 - 1d0);
        nvec(4) = l3*(2d0*l3 - 1d0);
        nvec(5) = 4d0 * l0 * l1;
        nvec(6) = 4d0 * l1 * l2;
        nvec(7) = 4d0 * l2 * l0;
        nvec(8) = 4d0 * l0 * l3;
        nvec(9) = 4d0 * l1 * l3;
        nvec(10) = 4d0 * l2 * l3;

        ftmp = 0.0d0
        do inode = 1, 10
            ftmp(3*inode-2) = nvec(inode) * load_direction(1)
            ftmp(3*inode-1) = nvec(inode) * load_direction(2)
            ftmp(3*inode) = nvec(inode) * load_direction(3)
        end do

        do inode = 1, 10
            rv(3*node_id(inode) - 2) = rv(3*node_id(inode) - 2) + ftmp(3*inode - 2)
            rv(3*node_id(inode) - 1) = rv(3*node_id(inode) - 1) + ftmp(3*inode - 1)
            rv(3*node_id(inode)) = rv(3*node_id(inode)) + ftmp(3*inode)
        end do
    end if

    call mpi_allreduce(found, found_cnt, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
    if (myid == 0) then
        print *, "number of process that found load element: ", found_cnt
    end if
    rv = rv / found_cnt

    ! output load_elem_arr
    write(filename, '("mdata/", I6.6, ".load_elem.bin")') myid
    open(10, file=filename, status="replace", access="stream", form="unformatted")
    write(10) load_elem_arr
    close(10)

end subroutine pointload

subroutine output_displacement(myid, n, ni, up)
    implicit none
    
    integer, intent(in) :: myid, n, ni
    double precision, intent(in) :: up(3*(n + ni))
    character(len=100) :: filename

    write(filename, '("displacement/", I6.6, ".bin")') myid
    open(10, file=filename, status="replace", access="stream", form="unformatted")
    write(10) up(1:3*n)
    close(10)

end subroutine output_displacement

subroutine check_load_elem(rv)
    implicit none

    double precision, intent(in) :: rv(:)
    integer :: i

    do i = 1, size(rv)
        if (rv(i) /= 0d0) then
            print *, "load element found!!"
            exit
        end if
    end do

end subroutine check_load_elem
