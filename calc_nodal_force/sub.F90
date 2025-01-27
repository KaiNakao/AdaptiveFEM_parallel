subroutine read_nload(nload)
    implicit none
    integer, intent(inout) :: nload
    character(len=80) :: filename

    write(filename, '(A)') 'data/pointload.dat'
    open(10, file=filename, status="old")
    read(10, *) ! number of loads
    read(10, *) nload
    close(10)
end subroutine read_nload


subroutine read_nobs_xnode(nobs, xnode)
    implicit none
    integer, intent(inout) :: nobs
    double precision, intent(inout) :: xnode(3, 8)
    character(len=80) :: filename
    integer :: iobs

    write(filename, '(A)') 'data/obs_points.dat'
    open(10, file=filename, status="old")
    read(10, *) ! number of observation points
    read(10, *) nobs
    read(10, *) ! x, y, z
    do iobs = 1, nobs - 8
        read(10, *) 
    end do
    print *, "xnode"
    do iobs = 1, 8
        read(10, *) xnode(:, iobs)
        print *, xnode(:, iobs)
    end do
    close(10)
end subroutine read_nobs_xnode

subroutine read_displacement(uvec, iload, nobs)
    implicit none
    double precision, intent(inout) :: uvec(24)
    integer, intent(in) :: iload, nobs
    double precision :: buf(3, nobs)
    integer :: inode
    character(len=80) :: filename        

    write(filename, '(A, I4.4, A)') 'displacement/', iload, '_obs.bin'
    open(10, file=filename, status="old", access="stream", form="unformatted")
    read(10) buf
    close(10)
    
    print *, "size(buf) = ", size(buf)

    do inode = 1, 8
        uvec(3*(inode - 1) + 1) = buf(1, nobs - 8 + inode)
        uvec(3*(inode - 1) + 2) = buf(2, nobs - 8 + inode)
        uvec(3*(inode - 1) + 3) = buf(3, nobs - 8 + inode)
    end do
    print *, "uvec = ", uvec
end subroutine read_displacement

subroutine read_centroid(xc, mvec)
    implicit none
    double precision, intent(inout) :: xc(3), mvec(6)
    character(len=80) :: filename

    write(filename, '(A)') 'data/target_centroid.dat'
    open(10, file=filename, status="old")
    read(10, *) ! centroid coordinate
    read(10, *) xc(1)
    read(10, *) xc(2)
    read(10, *) xc(3)
    read(10, *) ! moment tensor
    read(10, *) mvec(1)
    read(10, *) mvec(2)
    read(10, *) mvec(3)
    read(10, *) mvec(4)
    read(10, *) mvec(5)
    read(10, *) mvec(6)
    close(10)
end subroutine

subroutine calc_r(xnode, xc, r)
    implicit none
    double precision, intent(in) :: xnode(3, 8), xc(3)
    double precision, intent(out) :: r(3)    
    double precision :: x1, x2, y1, y2, z1, z2
    integer :: inode

    x1 = xnode(1, 1)
    x2 = xnode(1, 2)
    y1 = xnode(2, 1)
    y2 = xnode(2, 4)
    z1 = xnode(3, 1)
    z2 = xnode(3, 5)

    r(1) = (2d0 * xc(1) - x1 - x2) / (x2 - x1)
    r(2) = (2d0 * xc(2) - y1 - y2) / (y2 - y1)
    r(3) = (2d0 * xc(3) - z1 - z2) / (z2 - z1)
end subroutine calc_r

subroutine calc_fvec(fvec, xnode, xc, mvec)
    implicit none
    double precision, intent(inout) :: fvec(24)
    double precision, intent(in) :: xnode(3, 8), xc(3), mvec(6)
    double precision :: bmat(6, 24), r(3)

    call calc_r(xnode, xc, r) 
    call calc_bmat(bmat, r, xnode)
    fvec = matmul(transpose(bmat), mvec)

end subroutine calc_fvec

subroutine calc_bmat(bmat, r, xnode)
    implicit none
    double precision, intent(in) :: r(3), xnode(3, 8)
    double precision, intent(out) :: bmat(6, 24)
    integer :: inode
    double precision :: dndr(3, 8), dndx(3, 8), &
                        x1, x2, y1, y2, z1, z2, r1, r2, r3
    
    r1 = r(1)
    r2 = r(2)
    r3 = r(3)

    dndr(1, 1) = -0.125d0*(1.0d0 - r2)*(1.0d0 - r3)
    dndr(1, 2) = 0.125d0*(1.0d0 - r2)*(1.0d0 - r3)
    dndr(1, 3) = 0.125d0*(1.0d0 + r2)*(1.0d0 - r3)
    dndr(1, 4) = -0.125d0*(1.0d0 + r2)*(1.0d0 - r3)
    dndr(1, 5) = -0.125d0*(1.0d0 - r2)*(1.0d0 + r3)
    dndr(1, 6) = 0.125d0*(1.0d0 - r2)*(1.0d0 + r3)
    dndr(1, 7) = 0.125d0*(1.0d0 + r2)*(1.0d0 + r3)
    dndr(1, 8) = -0.125d0*(1.0d0 + r2)*(1.0d0 + r3)

    dndr(2, 1) = -0.125d0*(1.0d0 - r1)*(1.0d0 - r3)
    dndr(2, 2) = -0.125d0*(1.0d0 + r1)*(1.0d0 - r3)
    dndr(2, 3) = 0.125d0*(1.0d0 + r1)*(1.0d0 - r3)
    dndr(2, 4) = 0.125d0*(1.0d0 - r1)*(1.0d0 - r3)
    dndr(2, 5) = -0.125d0*(1.0d0 - r1)*(1.0d0 + r3)
    dndr(2, 6) = -0.125d0*(1.0d0 + r1)*(1.0d0 + r3)
    dndr(2, 7) = 0.125d0*(1.0d0 + r1)*(1.0d0 + r3)
    dndr(2, 8) = 0.125d0*(1.0d0 - r1)*(1.0d0 + r3)

    dndr(3, 1) = -0.125d0*(1.0d0 - r1)*(1.0d0 - r2)
    dndr(3, 2) = -0.125d0*(1.0d0 + r1)*(1.0d0 - r2)
    dndr(3, 3) = -0.125d0*(1.0d0 + r1)*(1.0d0 + r2)
    dndr(3, 4) = -0.125d0*(1.0d0 - r1)*(1.0d0 + r2)
    dndr(3, 5) = 0.125d0*(1.0d0 - r1)*(1.0d0 - r2)
    dndr(3, 6) = 0.125d0*(1.0d0 + r1)*(1.0d0 - r2)
    dndr(3, 7) = 0.125d0*(1.0d0 + r1)*(1.0d0 + r2)
    dndr(3, 8) = 0.125d0*(1.0d0 - r1)*(1.0d0 + r2)

    x1 = xnode(1, 1)
    x2 = xnode(1, 2)
    y1 = xnode(2, 1)
    y2 = xnode(2, 4)
    z1 = xnode(3, 1)
    z2 = xnode(3, 5)
    
    do inode = 1, 8
        dndx(1, inode) = dndr(1, inode) * 2d0 / (x2 - x1)
        dndx(2, inode) = dndr(2, inode) * 2d0 / (y2 - y1)
        dndx(3, inode) = dndr(3, inode) * 2d0 / (z2 - z1)
    end do

    bmat = 0d0
    do inode = 1, 8
        bmat(1, 3*(inode - 1) + 1) = dndx(1, inode)
        bmat(2, 3*(inode - 1) + 2) = dndx(2, inode)
        bmat(3, 3*(inode - 1) + 3) = dndx(3, inode)
        bmat(4, 3*(inode - 1) + 1) = dndx(2, inode)
        bmat(4, 3*(inode - 1) + 2) = dndx(1, inode)
        bmat(5, 3*(inode - 1) + 2) = dndx(3, inode)
        bmat(5, 3*(inode - 1) + 3) = dndx(2, inode)
        bmat(6, 3*(inode - 1) + 1) = dndx(3, inode)
        bmat(6, 3*(inode - 1) + 3) = dndx(1, inode)
    end do
end subroutine calc_bmat    

subroutine output_gf(gf)
    implicit none
    double precision, intent(in) :: gf(:)
    integer :: iload
    character(len=80) :: filename

    write(filename, '(A)') 'data/greens_function.dat'
    open(10, file=filename, status="replace")
    do iload = 1, size(gf)
        write(10, *) gf(iload)
    end do
    close(10)
end subroutine output_gf
