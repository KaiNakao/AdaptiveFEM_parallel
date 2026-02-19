#include "pointload_solver.hpp"

Pointload_Solver::Pointload_Solver(int myid, int numprocs) {
    // constructor
    im = myid;
    ncpu = numprocs;
}

void Pointload_Solver::prep() {
    initialize_(&im, &ncpu, &np, &kd, &nfa, &n2dvisl, &n2dvisa, 
                &n_size, &ne_size, &n, &ne, &ib, &nab, &nad, &nadto, 
                &n4_size, &ne4_size, &n4, &ne4, &ib4, &nab4, &nad4, &nadto4, 
                &bufls, &bufds, &fidc, &fidm, &ierr);

    plane.resize(5);
    coor.resize(3*n_size);
    // coor4.resize(3*n4_size);

    main_compute_ni_(&im, &ncpu, &np, &kd, &nfa, &n2dvisl,
            &n_size, &ne_size, &n, &ne, &ib, &nab, &nad, &nadto, &fidm,
            &n4_size, &ne4_size, &n4, &ne4, &ib4, &nab4, &nad4, &nadto4, &fidc,
            &bufds, &bufls,
            plane.data(), &ds, &ni_size, &nei_size, &ni, &nei, &ni4_size, &nei4_size,
            &ni4, &nei4);
    
    plane.resize(5);
    coor.resize(3*n_size); coori.resize(3*(n_size+ni_size));
    uval.resize(6*(n_size+ni_size));
    bcarray.resize(n_size+ni_size);
    cny.resize(10*ne_size); num.resize(ne_size);
    cnyi.resize(12*nei_size); numi.resize(nei_size); nfla.resize(nei_size);
    // dupli.resize(n_size); 
    duplixyz.resize(3*(n_size+ni_size));
    // ibi.resize(ib); 
    // nabi.resize(6*nab); 
    // npl.resize(np+1); 
    // abcelem.resize(nab);
    // mpiad.resize(nad*2); 
    // mpiadlist.resize(nadto); 
    // mpl.resize(nad+1);
    mpiadi.resize(nad*2); mpiadlisti.resize(2*nadto); mpli.resize(nad+1);
    // coors.resize(3*n_size);

    // coor4.resize(3*n4_size); 
    coori4.resize(3*(n4_size+ni4_size));
    uval4.resize(6*(n4_size+ni4_size));
    bcarray4.resize(n4_size+ni4_size);
    cny4.resize(4*ne4_size); num4.resize(ne4_size);
    cnyi4.resize(6*nei4_size); numi4.resize(nei4_size); nfla4.resize(nei4_size);
    // dupli4.resize(n4_size); 
    duplixyz4.resize(3*(n4_size+ni4_size));
    // ibi4.resize(ib4); nabi4.resize(3*nab4); npl4.resize(np+1); abcelem4.resize(nab4);
    // mpiad4.resize(nad4*2); mpiadlist4.resize(nadto4); mpl4.resize(nad4+1);
    mpiadi4.resize(nad4*2); mpiadlisti4.resize(2*nadto4); 
    mpli4.resize(nad4+1);

    pentawptrd.resize(2*nsta+3); pentarptrd.resize((n+ni)*4);
    pentawindd.resize((n+ni)*4); pentarindd.resize((n+ni)*CRSSIZE);
    pentawptr.resize(2*nsta+3); pentarptr.resize((n+ni)*4);
    pentawind.resize((n+ni)*4); pentarind.resize((n+ni)*CRSSIZE);
    pentawptr4.resize(2*nsta+3); pentarptr4.resize((n4+ni4)*4);
    pentawind4.resize((n4+ni4)*4); pentarind4.resize((n4+ni4)*CRSSIZE);

#if npc == 1 || npc == 21
    kcrsvald.resize(9*(n+ni)*CRSSIZE);
    kcrsval.resize(9*(n+ni)*CRSSIZE);
    kcrsval4.resize(9*(n4+ni4)*CRSSIZE);
#elif npc == 288
    kcrsvald.resize(9*(n+ni)*CRSSIZE*2);
    kcrsval.resize(9*(n+ni)*CRSSIZE*2);
    kcrsval4.resize(9*(n4+ni4)*CRSSIZE*2);
#endif
    map4to10.resize(2*(n+ni));
    
    // cnysep_read.resize(10*ne_size); 
    cnysep_write.resize(10*ne_size);
    // cnysep_read4.resize(4*ne4_size); 
    cnysep_write4.resize(4*ne4_size);
    nsep_ne_ptr.resize(nsta+2); numsep.resize(ne_size);
    nsep_n_ptr.resize(nsta+2);
    nsep_n_ptr4.resize(nsta+2);
    nsepnmap.resize(2*n);
    nsepnmap4.resize(2*n4);
    // numsep4.resize(ne); elemmap.resize(ne);
    color_ind.resize(np*MAXCOLOR);
    cnysep_readm.resize(10*ne); cnysep_readm4.resize(4*ne);

    // rmat.resize(10*kd); 
    younglst.resize(2*kd);
    // cort.resize(3); rf.resize(3); ftmp.resize(30); xx.resize(4*3);
    // corf.resize(3*nfa); faultp.resize(4*nfa); 
    sig.resize(nkl);
    up.resize(3*(n+ni)*npc); rv.resize(3*(n+ni)*npc);
    // u2dl.resize(3*n2dvisl*npc); utmp.resize(3*(n+ni)); utmp2.resize(3*(n+ni));
    // me.resize(10,10);
    // cori.resize(3*3); keit.resize(36*36); keit4.resize(18*18);
    // flag.resize(n_size);
    // u2dflagl.resize(n2dvisl);
    mpibuf.resize(2*3*npc*(nadto+1)*2);
    // node_old2new.resize(n); node_old2new4.resize(n4);
    uvalcoe.resize(npc);

    main_solver_(&im, &ncpu, &np, &kd, &nfa, plane.data(), &ds, &n2dvisl,
                &n_size, &ne_size, &ni_size, &nei_size, &n, &ne, &ni, &nei,
                &ib, &nab, &nad, &nadto, &fidm,
                &n4_size, &ne4_size, &ni4_size, &nei4_size, &n4, &ne4, &ni4, &nei4,
                &ib4, &nab4, &nad4, &nadto4, &fidc,
                &bufds, &bufls,
                coor.data(), coori.data(), &dsi, uval.data(), bcarray.data(), cny.data(), num.data(),
                cnyi.data(), numi.data(), nfla.data(), duplixyz.data(),
                mpiadi.data(), mpiadlisti.data(), mpli.data(), &nadtoi,
                coori4.data(), &dsi4, uval4.data(), bcarray4.data(), cny4.data(), num4.data(),
                cnyi4.data(), numi4.data(), nfla4.data(), duplixyz4.data(),
                mpiadi4.data(), mpiadlisti4.data(), mpli4.data(), &nadtoi4,
                &nird, &niwd, pentawptrd.data(), pentarptrd.data(), pentawindd.data(), pentarindd.data(),
                &nir, &niw, pentawptr.data(), pentarptr.data(), pentawind.data(), pentarind.data(),
                &nir4, &niw4, pentawptr4.data(), pentarptr4.data(), pentawind4.data(), pentarind4.data(),
                kcrsvald.data(), kcrsval.data(), kcrsval4.data(),
                map4to10.data(),
                cnysep_write.data(), cnysep_write4.data(),
                nsep_ne_ptr.data(), numsep.data(), nsep_n_ptr.data(), nsep_n_ptr4.data(), &nsep,
                nsepnmap.data(), &nsep4, nsepnmap4.data(),
                &ncolor, color_ind.data(),
                cnysep_readm.data(), cnysep_readm4.data(),
                younglst.data(), sig.data(), up.data(), rv.data(), mpibuf.data(), uvalcoe.data(),
                &log_innerCG, &log_innerCG_per_outite,
                &use_crs_d, &overlap_comm_d,
                &overlap_comm_s, &persistent_comm_s);

    if (im == 0) {
        nload = 0;
        read_nload_(&nload);
    }
    MPI_Bcast(&nload, 1, MPI_INT, 0, MPI_COMM_WORLD);
    load_arr.resize(5*nload);
    if (im == 0) {
        read_loads_(&nload, load_arr.data());
    }
    MPI_Bcast(load_arr.data(), 5*nload, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    surf_nelem = 0;
    read_surf_nelem_(&im, &surf_nelem);
    if (surf_nelem > 0) {
        surf_cny.resize(4*surf_nelem);
        read_surf_cny_(&im, &surf_nelem, surf_cny.data());
    }

    load_elem_arr.resize(ne, 0);
}

void Pointload_Solver::exec(int iload, bool write_displacement) {
    int iload_fortran = iload + 1; // Fortran is 1-indexed
    pointload_(&im, &n, coor.data(), &ne, cny.data(), &ni, rv.data(), 
              load_arr.data(), &nload, &iload_fortran, load_elem_arr.data(), 
              &surf_nelem, surf_cny.data());
    int nni = n + ni;
    sync_parallel_block_d_(&nni,&nni, 
            rv.data(),&nad,&nadtoi,mpiadi.data(),mpiadlisti.data(),mpli.data(),&ierr,&im,mpibuf.data());
    for (int i = 0; i < up.size(); ++i) {
        up.at(i) = 0.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double t1=MPI_Wtime();
    ipcg_(&im,&np,&kd,younglst.data(),sig.data(),mpibuf.data(),
          &n,&ne,coori.data(),cny.data(),num.data(),uval.data(),duplixyz.data(),bcarray.data(),
          &ni,&nei,cnyi.data(),&dsi,numi.data(),nfla.data(),
          &nad,&nadtoi,mpiadi.data(),mpiadlisti.data(),mpli.data(),
          &nird,pentarptrd.data(),pentarindd.data(),
          &niwd,pentawptrd.data(),pentawindd.data(),kcrsvald.data(),
          &nir,pentarptr.data(),pentarind.data(),
          &niw,pentawptr.data(),pentawind.data(),kcrsval.data(),
          &nsep,numsep.data(),nsep_n_ptr.data(),nsep_ne_ptr.data(),nsepnmap.data(),
          cnysep_readm.data(),cnysep_write.data(),
          &ncolor,color_ind.data(),
          &n4,&ne4,coori4.data(),cny4.data(),num4.data(),uval4.data(),duplixyz4.data(),bcarray4.data(),
          &ni4,&nei4,cnyi4.data(),&dsi4,numi4.data(),nfla4.data(),
          &nad4,&nadtoi4,mpiadi4.data(),mpiadlisti4.data(),mpli4.data(),
          &nir4,pentarptr4.data(),pentarind4.data(),
          &niw4,pentawptr4.data(),pentawind4.data(),kcrsval4.data(),
          map4to10.data(),
          cnysep_readm4.data(),cnysep_write4.data(),nsep_n_ptr4.data(),&nsep4,nsepnmap4.data(),
          uvalcoe.data(),&log_innerCG,&log_innerCG_per_outite,
          &use_crs_d,&overlap_comm_d,&overlap_comm_s,&persistent_comm_s,
          rv.data(),up.data());
    MPI_Barrier(MPI_COMM_WORLD);
    double t2=MPI_Wtime();
    if (im == 0) {
        std::cout << "IPCG took " << t2 - t1 << " seconds." << std::endl;
    }
    if (write_displacement) {
        output_displacement_(&im, &n, &ni, up.data(), &iload_fortran);
    }
}
