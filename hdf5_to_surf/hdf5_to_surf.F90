! -s オプションはモデルに空洞がないかを確認するためのもの
module domain_model
  implicit none
  type FEmodel
    character(len = 100) :: name
    integer(8) :: n, ne, kd, npe = 4
    integer(8), allocatable :: cny(:, :)
    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(8), allocatable :: coor(:, :)
  end type FEmodel
contains
subroutine allocate_model(model)
  implicit none
  type(FEmodel), intent(inout) :: model
  allocate(model%cny(model%npe+1, model%ne))
  allocate(model%coor(3, model%n))
end subroutine allocate_model

subroutine deallocate_model(model)
  implicit none
  type(FEmodel), intent(inout) :: model
  deallocate(model%cny, model%coor)
end subroutine deallocate_model

subroutine calc_domain_range(model)
  implicit none
  type(FEmodel), intent(inout) :: model
  integer(8) :: i
  real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
  if(.not. allocated(model%coor)) stop "The model is not allocated in calc_domain_range."
  if(model%n == 0) return
  i = 1
  xmin = model%coor(1,i)
  xmax = model%coor(1,i)
  ymin = model%coor(2,i)
  ymax = model%coor(2,i)
  zmin = model%coor(3,i)
  zmax = model%coor(3,i)
  do i = 1, model%n
    xmin = min(xmin, model%coor(1,i))
    xmax = max(xmax, model%coor(1,i))
    ymin = min(ymin, model%coor(2,i))
    ymax = max(ymax, model%coor(2,i))
    zmin = min(zmin, model%coor(3,i))
    zmax = max(zmax, model%coor(3,i))
  enddo
  write(6, *) "domain range of model ", trim(model%name)
  write(6, *) "--- x", xmin, xmax
  write(6, *) "--- y", ymin, ymax
  write(6, *) "--- z", zmin, zmax
  model%xmin = xmin
  model%xmax = xmax
  model%ymin = ymin
  model%ymax = ymax
  model%zmin = zmin
  model%zmax = zmax
end subroutine calc_domain_range

subroutine print_aspect_ratio(model)
  implicit none
  type(FEmodel), intent(inout) :: model
  integer(8) :: ie, i1, i2, i3, i4
  real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3, vol, vol0, radius, ratio
  real(8) :: a, b, c, d, e, f, ad, be, cf
  real(8) :: rmax, rmin, Vsum, Vmin
  integer(8) :: ne_negative_vol, ne_small_aratio
  rmax = -1d100
  rmin = 1d100
  Vsum = 0.d0
  Vmin = 1d100
  ne_negative_vol = 0
  ne_small_aratio = 0
  do ie = 1, model%ne
    !swap = 0
    i1 = model%cny(1, ie)
    i2 = model%cny(2, ie)
    i3 = model%cny(3, ie)
    i4 = model%cny(4, ie)
    x1 = model%coor(1, i2) - model%coor(1, i1)
    x2 = model%coor(1, i3) - model%coor(1, i1)
    x3 = model%coor(1, i4) - model%coor(1, i1)
    y1 = model%coor(2, i2) - model%coor(2, i1)
    y2 = model%coor(2, i3) - model%coor(2, i1)
    y3 = model%coor(2, i4) - model%coor(2, i1)
    z1 = model%coor(3, i2) - model%coor(3, i1)
    z2 = model%coor(3, i3) - model%coor(3, i1)
    z3 = model%coor(3, i4) - model%coor(3, i1)
    vol = x1*y2*z3+x2*y3*z1+x3*y1*z2-x1*y3*z2-x2*y1*z3-x3*y2*z1
    if(vol <= 0.)then
      ne_negative_vol = ne_negative_vol + 1
      model%cny(5,ie) = -1
    endif
    vol = vol/6.d0
    Vsum = Vsum + vol
    Vmin = min(Vmin,vol)
    a = sqrt(x1**2+y1**2+z1**2)
    b = sqrt(x2**2+y2**2+z2**2)
    c = sqrt(x3**2+y3**2+z3**2)
    d = sqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
    e = sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    f = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    ad = a*d
    be = b*e
    cf = c*f
    !外接球の半径
    if((ad+be+cf)*(-ad+be+cf)*(ad-be+cf)*(ad+be-cf)<0. .or. abs(vol)<1e-20)then
      write(6,*)ie,(ad+be+cf)*(-ad+be+cf)*(ad-be+cf)*(ad+be-cf),vol
    endif
    radius = sqrt((ad+be+cf)*(-ad+be+cf)*(ad-be+cf)*(ad+be-cf))/(24*vol)
    !外接球の半径がrであるような正四面体の体積
    vol0 = 8/(9*sqrt(3.))*radius**3
    ! アスペクト比
    ratio = vol/vol0
    if(ratio<0.006)then
      ne_small_aratio = ne_small_aratio
      model%cny(5,ie) = 0
    endif
    rmin = min(rmin, ratio)
    rmax = max(rmax, ratio)
  enddo
  write(6, *) "Aspect ratio of model ", trim(model%name)
  if(ne_small_aratio > 0) write(6, *) "WARNING: num of elem. with too small aspect ratio", ne_small_aratio
  if(ne_negative_vol > 0) write(6, *) "WARNING: num of elem. with negative volume", ne_negative_vol
  write(6,*)"--- range of aspect ratio:",rmin,rmax
  write(6,*)"--- Volume               :",Vsum,Vmin
end subroutine print_aspect_ratio

subroutine read_model_from_hdf(model)
  implicit none
  type(FEmodel), intent(inout) :: model
  character(len = 100) :: filename
  integer(8) :: n, ne, ne_vox, kd
  integer(8) :: ltmp, fid, bufsize, datasize
  !------------------------------------------------------------------------------------
  ! read setting
  filename = "./data/tet4vox8setting" // trim(model%name) // ".dat"
  open(10, file = filename, status = 'old')
  read(10, *)
  read(10, *) n
  read(10, *)
  read(10, *) ne
  read(10, *)
  read(10, *) ne_vox
  read(10, *)
  read(10, *) kd
  close(10)
  model%n = n
  model%ne = ne
  model%kd = kd
  if(ne_vox > 0) stop 'This code only supports ne_vox==0'
  write(6,*) "n     :",n
  write(6,*) "ne    :",ne
  write(6,*) "ne_vox:",ne_vox
  write(6,*) "kd    :",kd
  !------------------------------------------------------------------------------------
  call allocate_model(model)
  !------------------------------------------------------------------------------------
  ! read hdf5 data
  fid = 0
  filename = "./data/tet4vox8model" // trim(model%name) // ".h5"
  write(6, *) 'reading model from ', trim(filename)
  ltmp = len(trim(filename))
  call set_null_char(filename, ltmp)
  call hdf_open_file(filename, fid)
  ! cny
  ltmp=10
  bufsize = 5 * ne
  write(filename, '(a10)') "/conn_tet4"
  call set_null_char(filename, ltmp)
  call hdf_read_long_array(fid, filename, model%cny, bufsize, datasize)
  if(datasize /= ne*5) stop "datasize differs : /conn_tet4"
  ! coor
  ltmp=5
  bufsize = 3 * n
  write(filename, '(a5)')"/coor"
  call set_null_char(filename, ltmp)
  call hdf_read_double_array(fid, filename, model%coor, bufsize, datasize)
  if(datasize /= n*3) stop "datasize differs : /coor"
  call hdf_close_file(fid)
  !end read hdf5 data
  !------------------------------------------------------------------------------------
end subroutine read_model_from_hdf

subroutine write_model_to_hdf(model)
  implicit none
  type(FEmodel), intent(in) :: model
  character(len = 100) :: filename
  integer(8) :: ltmp, fid, bufsize
  !------------------------------------------------------------------------------------
  ! write setting
  filename = "./data/tet4vox8setting" // trim(model%name) // ".dat"
  open(10, file = filename, status = 'unknown')
  write(10, *) 'num of node'
  write(10, *) model%n
  write(10, *) 'num of tet elem'
  write(10, *) model%ne
  write(10, *) 'num of vox elem'
  write(10, *) 0
  write(10, *) 'num of material properties'
  write(10, *) model%kd
  close(10)
  !------------------------------------------------------------------------------------
  ! read hdf5 data
  fid = 0
  filename = "./data/tet4vox8model" // trim(model%name) // ".h5"
  write(6, *) 'writing model to ', trim(filename)
  ltmp = len(trim(filename))
  call set_null_char(filename, ltmp)
  call hdf_create_file(filename, fid)
  ! cny
  ltmp=10
  bufsize = 5 * model%ne
  write(filename, '(a10)') "/conn_tet4"
  call set_null_char(filename, ltmp)
  call hdf_write_long_array(fid, filename, model%cny, bufsize, bufsize)
  ! coor
  ltmp=5
  bufsize = 3 * model%n
  write(filename, '(a5)') "/coor"
  call set_null_char(filename, ltmp)
  call hdf_write_double_array(fid, filename, model%coor, bufsize, bufsize)
  call hdf_close_file(fid)
  !end read hdf5 data
  !------------------------------------------------------------------------------------
end subroutine write_model_to_hdf

end module domain_model

module sort_array
contains
subroutine sort(x)
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(8), intent(inout) :: x(:, :)
  integer(8) :: xsize
  integer(8), allocatable :: pivot(:), tmp(:)
  xsize = ubound(x, 1)
  allocate(pivot(xsize), tmp(xsize))
  call quick_sort(1, ubound(x, 2, int64))
  deallocate(pivot, tmp)
contains
  recursive subroutine quick_sort(ista, iend)
    implicit none
    integer(8), intent(in) :: ista, iend
    integer(8) :: ismall, ilarge
    if(iend <= ista) return
    pivot(:) = x(:, ista)
    ismall = ista
    ilarge = iend
    do while(ismall < ilarge)
      do while(.not. is_larger(x(:, ismall), pivot) .and. ismall <= iend)
        ismall = ismall + 1
      enddo
      do while(is_larger(x(:, ilarge), pivot) .and. ilarge >= ista)
        ilarge = ilarge - 1
      enddo
      if(ilarge <= ismall) exit
      tmp = x(:, ismall)
      x(:, ismall) = x(:, ilarge)
      x(:, ilarge) = tmp
    enddo
    ismall = ismall - 1
    x(:, ista) = x(:, ismall)
    x(:, ismall) = pivot
    call quick_sort(ista, ismall - 1)
    call quick_sort(ismall + 1, iend)
  end subroutine quick_sort
end subroutine sort

logical function is_larger(x1, x2)
  implicit none
  integer(8), intent(in) :: x1(:), x2(:)
  integer(8) :: xsize, j1
  xsize = ubound(x1, 1)
  is_larger = .false.
  do j1 = 1, xsize
    if(x1(j1) > x2(j1))then
      is_larger = .true.
      exit
    elseif(x1(j1) < x2(j1))then
      is_larger = .false.
      exit
    endif
  enddo
end function is_larger
end module sort_array

program main
  use domain_model
  implicit none
  type(FEmodel) :: model, model_surf
  logical :: no_material_boundary, add_split_id

  call setting(model, model_surf, no_material_boundary, add_split_id)
  call read_model_from_hdf(model)
  call extract_surface_elements(model, model_surf, no_material_boundary, add_split_id)
  call write_model_to_hdf(model_surf)

contains
subroutine setting(model, model_surf, no_material_boundary, add_split_id)
  implicit none
  type(FEmodel), intent(out) :: model, model_surf
  logical, intent(out) :: no_material_boundary, add_split_id
  character(len = 50) :: ctmp
  integer(8) :: iarg_model_name = 0
  iarg_model_name = 1
  no_material_boundary = .false.
  add_split_id = .false.
  if(command_argument_count() > 0)then
    call getarg(1, ctmp)
    if(trim(ctmp) == 'help')then
      write(6, *) './hdf5_to_surf.exe (-n) (-s) ([model name])'
      write(6, *) '   -n: without material boundary.'
      write(6, *) '   -s: output surface as vtk for each split'
      write(6, *) '   "tet4vox8model[model name].h5" is the target file.'
      write(6, *) './hdf5_to_surf.exe help: show this help'
      stop
    elseif(trim(ctmp) == '-s')then
      no_material_boundary = .true.
      add_split_id = .true.
      iarg_model_name = 2
    elseif(trim(ctmp) == '-n')then
      no_material_boundary = .true.
      iarg_model_name = 2
    endif
  endif
  if(command_argument_count() < iarg_model_name)then
    model%name = ''
  else
    call getarg(iarg_model_name, model%name)
  endif

  model_surf%name = trim(model%name) // '.surf'
  write(6, *) "target model: ", trim(model%name)
  write(6, *) "no_material_boundary:", no_material_boundary
  write(6, *) "add_split_id:", add_split_id
end subroutine setting

subroutine extract_surface_elements(model, model_surf, no_material_boundary, add_split_id)
  use sort_array
  implicit none
  type(FEmodel), intent(inout) :: model, model_surf
  logical, intent(in) :: no_material_boundary, add_split_id
  integer(8) :: tri_surf(4, model%ne*4)
  integer(8) :: cny_sorted(4)
  integer(8) :: ie, j, i, ne_surf, n_surf
  integer(8) :: org2surf(model%n)
  logical :: elem_on_surf(model%ne), node_on_surf(model%n)
  integer(8) :: n_split, tri_surf_split_id(model%ne*4), elem_surf_split_id(model%ne)

  ! set tri_surf
  do ie = 1, model%ne
    j = (ie-1) * 4
    call sort_cny(model%cny(1:4, ie), cny_sorted)
    tri_surf(1:3, j + 1) = (/cny_sorted(1), cny_sorted(2), cny_sorted(3)/)
    tri_surf(1:3, j + 2) = (/cny_sorted(1), cny_sorted(2), cny_sorted(4)/)
    tri_surf(1:3, j + 3) = (/cny_sorted(1), cny_sorted(3), cny_sorted(4)/)
    tri_surf(1:3, j + 4) = (/cny_sorted(2), cny_sorted(3), cny_sorted(4)/)
    tri_surf(4, j+1:j+4) = ie
  enddo

  ! sort tri_surf
  call sort(tri_surf)

  ! remove non surface tri
  j = 1
  do while(j < model%ne * 4)
    if(all(tri_surf(1:3, j) == tri_surf(1:3, j + 1)))then ! 同じ3角形パッチ
      if(model%cny(5, tri_surf(4, j)) == model%cny(5, tri_surf(4, j + 1)) & ! 要素の物性が同じなら
      & .or. no_material_boundary)then ! または、 no_material_boundary == .true. なら
        ! set removed flag
        tri_surf(4, j) = -1
        tri_surf(4, j + 1) = -1
      endif
      j = j + 1
    endif
    j = j + 1
  enddo

  if(add_split_id)then
    call set_tri_split_id(model%n, model%ne*4, tri_surf, n_split, tri_surf_split_id)
    call reorder_tri_node(model, tri_surf)
    call output_tri_surf_split(model, n_split, tri_surf, tri_surf_split_id)
  endif

  elem_on_surf(:) = .false.
  do j = 1, model%ne * 4
    ie = tri_surf(4, j)
    if(ie > 0)then
      elem_on_surf(ie) = .true.
    endif
  enddo

  ne_surf = 0
  node_on_surf(:) = .false.
  do ie = 1, model%ne
    if(elem_on_surf(ie))then
      ne_surf = ne_surf + 1
      do j = 1, 4
        i = model%cny(j, ie)
        node_on_surf(i) = .true.
      enddo
    endif
  enddo
  n_surf = 0
  do i = 1, model%n
    if(node_on_surf(i))then
      n_surf = n_surf + 1
      org2surf(i) = n_surf
    endif
  enddo

  model_surf%n = n_surf
  model_surf%ne = ne_surf
  model_surf%kd = model%kd
  call allocate_model(model_surf)

  j = 0
  do i = 1, model%n
    if(node_on_surf(i))then
      j = j + 1
      model_surf%coor(:, j) = model%coor(:, i)
    endif
  enddo

  j = 0
  do ie = 1, model%ne
    if(elem_on_surf(ie))then
      j = j + 1
      model_surf%cny(1:4, j) = org2surf(model%cny(1:4, ie))
      model_surf%cny(5, j) = model%cny(5, ie)
    endif
  enddo
end subroutine extract_surface_elements

subroutine set_tri_split_id(n, ntri, tri, n_split, tri_split_id)
  implicit none
  integer(8), intent(in) :: n, ntri, tri(4, ntri)
  integer(8), intent(out) :: n_split, tri_split_id(ntri)
  integer(8) :: ntri_used, split_id_node(n)
  integer(8), allocatable :: tri_used(:,:), used2all(:), tri_used_split_id(:)
  integer(8), allocatable :: split_id_tri(:), split_id_old2new(:)
  integer(8) :: i, j1, itri, itri_used, id_tri, id_node, loop, count_changed
  integer, parameter :: NOT_SET = -1
  ! tri: 3角形。tri(4,i) は三角形iの所属する要素。
  ! ただし、この時点では、モデルの表面にない三角形i については tri(4,i)=-1
  ! このsubroutine では、tri(4,i) > 0 の三角形のみが興味の対象

  ! 1. 使う三角形をtri_used に抜き出す
  !    tri_used(4, :) に split_id を入れる
  ntri_used = 0
  do itri = 1, ntri
    if(tri(4, itri) > 0) ntri_used = ntri_used + 1
  enddo
  allocate(tri_used(3, ntri_used))
  allocate(used2all(ntri_used), tri_used_split_id(ntri_used))
  itri_used = 0
  do itri = 1, ntri
    if(tri(4, itri) > 0)then
      itri_used = itri_used + 1
      tri_used(1:3, itri_used) = tri(1:3, itri)
      used2all(itri_used) = itri
    endif
  enddo

  ! 2. split_id
  allocate(split_id_tri(ntri_used))
  ! 2-1. split_id の初期化
  do itri = 1, ntri_used
    split_id_tri(itri) = itri
  enddo
  split_id_node(:) = NOT_SET
  ! 2-2. split_id の更新
  do loop = 1, 10
    count_changed = 0
    do itri = 1, ntri_used
      ! 各三角形について、要素と節点のsplit_idのうち最も小さいものを
      ! 要素と節点のsplit_id にする
      id_tri = split_id_tri(itri)
      do j1 = 1, 3
        i = tri_used(j1, itri)
        id_node = split_id_node(i)
        if(id_node > 0)then
          id_tri = min(id_tri, split_id_tri(id_node))
        endif
      enddo
      if(split_id_tri(itri) /= id_tri) count_changed = count_changed + 1
      split_id_tri(itri) = id_tri
      do j1 = 1, 3
        i = tri_used(j1, itri)
        id_node = split_id_node(i)
        split_id_node(i) = id_tri
        if(id_node /= id_tri) count_changed = count_changed + 1
        if(id_node >0)then
          if(split_id_tri(id_node) /= id_tri) count_changed = count_changed + 1
          split_id_tri(id_node) = id_tri
        endif
      enddo
    enddo
    write(6, *) "split loop", loop, count_changed
    if(count_changed == 0)then
      write(6, *) "splitted successfully."
      exit
    endif
  enddo ! loop

  ! 3. split_id_tri に使われている split_id を抽出し、1から順番に番号を振りなおす
  allocate(split_id_old2new(ntri_used))
  split_id_old2new(:) = NOT_SET
  do itri_used = 1, ntri_used
    id_tri = split_id_tri(itri_used)
    split_id_old2new(id_tri) = 1
  enddo
  id_tri = 0
  do itri_used = 1, ntri_used
    if(split_id_old2new(itri_used) == NOT_SET) cycle
    id_tri = id_tri + 1
    split_id_old2new(itri_used) = id_tri
  enddo
  n_split = id_tri
  write(6, *) "num of splits:", n_split
  do itri_used = 1, ntri_used
    id_tri = split_id_tri(itri_used)
    split_id_tri(itri_used) = split_id_old2new(id_tri)
  enddo

  ! 3. tri_split_id に書き出す
  tri_split_id(:) = NOT_SET
  do itri_used = 1, ntri_used
    itri = used2all(itri_used)
    tri_split_id(itri) = split_id_tri(itri_used)
  enddo
end subroutine set_tri_split_id

subroutine reorder_tri_node(model, tri_surf)
  ! tri_surf の節点の順番は、節点番号が昇順になるようになっているが、
  ! 今後使うときにはもとの順番になっていてほしいので、並び替える
  implicit none
  type(FEmodel), intent(in) :: model
  integer(8), intent(inout) :: tri_surf(4, model%ne*4)
  integer(8) :: itri, j1, j2, flag, i1, i2, ie, iflag
  do itri = 1, model%ne * 4
    if(tri_surf(4, itri) < 0) cycle
    ie = tri_surf(4, itri)
    flag = 1 + 2 + 4 + 8
    do j1 = 1, 4
      i1 = model%cny(j1, ie)
      do j2 = 1, 3
        i2 = tri_surf(j2, itri)
        if(i1 == i2)then
          flag = flag - 2 ** (j1-1)
          exit
        endif
      enddo ! j1
    enddo

    select case(flag)
    case(1)
      tri_surf(1, itri) = model%cny(2, ie)
      tri_surf(2, itri) = model%cny(4, ie)
      tri_surf(3, itri) = model%cny(3, ie)
    case(2)
      tri_surf(1, itri) = model%cny(1, ie)
      tri_surf(2, itri) = model%cny(3, ie)
      tri_surf(3, itri) = model%cny(4, ie)
    case(4)
      tri_surf(1, itri) = model%cny(1, ie)
      tri_surf(2, itri) = model%cny(4, ie)
      tri_surf(3, itri) = model%cny(2, ie)
    case(8)
      tri_surf(1, itri) = model%cny(1, ie)
      tri_surf(2, itri) = model%cny(2, ie)
      tri_surf(3, itri) = model%cny(3, ie)
    case default
      write(6, *) flag
      write(6, *) tri_surf(1:3, itri)
      write(6, *) model%cny(1:4, ie)
      stop "ERROR: in reorder_tri_node."
    end select
  enddo ! itri
end subroutine reorder_tri_node

subroutine output_tri_surf_split(model, n_split, tri_surf, tri_surf_split_id)
  implicit none
  type(FEmodel), intent(in) :: model
  integer(8), intent(in) :: n_split, tri_surf(4, model%ne*4), tri_surf_split_id(model%ne*4)
  integer(8) :: n, ntri, node_split(model%n), tri_split(model%ne*4), all2split(model%n)
  integer(8) :: itri, i_split, j1, i, conn(3)
  character(len = 100) :: filename
  integer, parameter :: NOT_SET = -1

  do i_split = 1, n_split
    all2split(:) = NOT_SET
    ntri = 0
    do itri = 1, model%ne * 4
      if(tri_surf_split_id(itri) /= i_split) cycle

      ntri = ntri + 1
      tri_split(ntri) = itri
      do j1 = 1, 3
        i = tri_surf(j1, itri)
        all2split(i) = 1
      enddo
    enddo
    n = 0
    do i = 1, model%n
      if(all2split(i) == NOT_SET) cycle
      n = n + 1
      node_split(n) = i
      all2split(i) = n
    enddo

    write(filename, '(a,i3.3,a)') "data/tri_surf.", i_split, ".vtk"
    open(20, file = filename, status = 'unknown')
    write(20, '(a)') "# vtk DataFile Version 4.0"
    write(20, '(a)') "file"
    write(20, '(a)') "ASCII"
    write(20, '(a)') "DATASET UNSTRUCTURED_GRID"
    write(20, '(a, i, a)') "POINTS", n, " double"
    do i = 1, n
      write(20, *) model%coor(:, node_split(i))
    enddo
    write(20, '(a, 2i)') "CELLS", ntri, ntri * 4
    do itri = 1, ntri
      do j1 = 1, 3
        conn(j1) = all2split(tri_surf(j1, tri_split(itri))) - 1
      enddo
      write(20, *) 3, conn
    enddo
    write(20, '(a, i)') "CELL_TYPES", ntri
    do itri = 1, ntri
      write(20, '(a)') '5'
    enddo
    write(20, '(a, i)') "CELL_DATA", ntri
    write(20, '(a)') "SCALARS material int 1"
    write(20, '(a)') "LOOKUP_TABLE default"
    do itri = 1, ntri
      write(20, '(i)') model%cny(5, tri_surf(4, tri_split(itri)))
    enddo
    close(20)
  enddo
end subroutine output_tri_surf_split

subroutine sort_cny(cny, cny_sorted)
  implicit none
  integer(8), intent(in) :: cny(4)
  integer(8), intent(out) :: cny_sorted(4)
  integer(8) :: j1, j2, itmp
  cny_sorted(1:4) = cny(1:4)
  do j1 = 1, 4
    do j2 = j1 + 1, 4
      if(cny_sorted(j1) > cny_sorted(j2))then
        itmp = cny_sorted(j1)
        cny_sorted(j1) = cny_sorted(j2)
        cny_sorted(j2) = itmp
      endif
    enddo
  enddo
end subroutine sort_cny
end program main
