#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

void convert_4to10_grid_cpp( long* n4_, long* ne_, long* coorsize_, double* coor, long* conn,
long* nl_, double* dsg_,long* nxg_, long* nyg_,long* nzg_,double* co,long* n10_,long* flag_)
{
  long n4,ne,nxg,nyg,nzg,n10,coorsize;
  double dsg;
  long i,j,ii,k,ix,iy,iz,k1,j1,i1,nm1,nm2,i2;
  long cnyele[10][2];
  double cri,coorlocal[10][3];
  long istart, iend, jstart, jend, kstart, kend;
  bool flag;
  
  n4=n4_[0];
  ne=ne_[0];
  nxg=nxg_[0];
  nyg=nyg_[0];
  nzg=nzg_[0];
  coorsize=coorsize_[0];
  dsg=dsg_[0];
  flag_[0]=1;
  n10=n4;
  cri=0.001;

  cnyele[4][0] = 0; //(5,1)=1
  cnyele[4][1] = 1; //(5,2)=2
  cnyele[5][0] = 1; //(6,1)=2
  cnyele[5][1] = 2; //(6,2)=3
  cnyele[6][0] = 0; //(7,1)=1
  cnyele[6][1] = 2; //(7,2)=3
  cnyele[7][0] = 0; //(8,1)=1
  cnyele[7][1] = 3; //(8,2)=4
  cnyele[8][0] = 1; //(9,1)=2
  cnyele[8][1] = 3; //(9,2)=4
  cnyele[9][0] = 3; //(10,1)=4
  cnyele[9][1] = 2; //(10,2)=3
 
//  std::cout << "n4,ne,nxg,nyg,nzg: " << n4 << " " << ne << " " << nxg << " " << nyg << " "<< nzg << "\n";
//  std::cout << "coorsize: " << coorsize << "\n";
//  while(nxg*nyg*nzg > 100000000)
  while(nxg*nyg*nzg > 500000000)
  {
//    std::cout << "make dsg larger " << ne << "\n";
//    std::cout.flush();
    nxg = nxg/2+1;
    nyg = nyg/2+1;
    nzg = nzg/2+1;
    dsg*=2.0;
  }
  std::vector<int> tmplist(0);
  std::vector<std::vector<int> > listn(nxg*nyg*nzg,tmplist);
      
  for( i = 0; i < ne; ++i )
  {
    for( j = 0; j < 4; ++j )
    {
      for( ii = 0; ii < 3; ++ii )
        coorlocal[j][ii] = coor[3*(conn[11*i+j]-1)+ii];
    }
    for( ii = 0; ii < 3; ++ii)
    {
      for( k = 4; k < 10; ++k )
        coorlocal[k][ii] = 0.5*(coorlocal[cnyele[k][0]][ii]+coorlocal[cnyele[k][1]][ii]);
    }
    for( k = 4; k < 10; ++k)
    {
    /*
      ix=int((coorlocal[k][0]-co[0])/dsg)+1;
      if(ix + 1 - ((coorlocal[k][0]-co[0])/dsg + 1) < 0.001)
      {
        istart = ix;
        iend   = ix + 1;
      }
      else if ((coorlocal[k][0]-co[0])/dsg + 1 - ix < 0.001)
      {
        istart = ix - 1;
        iend   = ix;
      }
      else
      {
        istart = ix;
        iend   = ix;
      }
      istart=ix-1;
      iend  =ix+1;

      iy=int((coorlocal[k][1]-co[1])/dsg)+1;
      if(iy + 1 - ((coorlocal[k][1]-co[1])/dsg + 1) < 0.001)
      {
        jstart = iy;
        jend   = iy + 1;
      }
      else if ((coorlocal[k][1]-co[1])/dsg + 1 - iy < 0.001)
      {
        jstart = iy - 1;
        jend   = iy;
      }
      else
      {
        jstart = iy;
        jend   = iy;
      }
      jstart=iy-1;
      jend  =iy+1;

      iz=int((coorlocal[k][2]-co[2])/dsg)+1;
      if(iz + 1 - ((coorlocal[k][2]-co[2])/dsg + 1) < 0.001)
      {
        kstart = iz;
        kend   = iz + 1;
      }
      else if ((coorlocal[k][2]-co[2])/dsg + 1 - iz < 0.001)
      {
        kstart = iz - 1;
        kend   = iz;
      }
      else{
        kstart = iz;
        kend   = iz;
      }
      kstart=iz-1;
      kend  =iz+1;
      */
     
      ix = round((coorlocal[k][0]-co[0])/dsg);
      iy = round((coorlocal[k][1]-co[1])/dsg);
      iz = round((coorlocal[k][2]-co[2])/dsg);
      istart = (ix-1 < 0) ? 0 : ix-1;
      jstart = (iy-1 < 0) ? 0 : iy-1;
      kstart = (iz-1 < 0) ? 0 : iz-1;
      iend = (ix+1 > nxg-1) ? nxg-1 : ix+1;
      jend = (iy+1 > nyg-1) ? nyg-1 : iy+1;
      kend = (iz+1 > nzg-1) ? nzg-1 : iz+1;
      
      flag=true;

      for( k1 = kstart; k1 <= kend && flag; ++k1 )
      {
      for( j1 = jstart; j1 <= jend && flag; ++j1 )
      {
      for( i1 = istart; i1 <= iend && flag; ++i1 )
      {
        nm1 = listn[nyg*nxg*k1+nxg*j1+i1].size();
        for( i2 = 0; i2 < nm1 && flag; ++i2 )
        {
          nm2 = listn[nyg*nxg*k1+nxg*j1+i1][i2];
          if( (fabs(coorlocal[k][0]-coor[3*(nm2-1)+0]) < cri) &&
              (fabs(coorlocal[k][1]-coor[3*(nm2-1)+1]) < cri) &&
              (fabs(coorlocal[k][2]-coor[3*(nm2-1)+2]) < cri) )
          {
            conn[11*i+k] = nm2;
            flag=false;
          }
        }
      }
      }
      }
      
      if(flag)
      {
        if( (n10+2) == coorsize )
        {
          std::cout << "coorsize too small " << coorsize << "\n";
          exit(0);
        }
        n10++;
        conn[11*i+k] = n10;
        for( ii = 0; ii < 3; ++ii )
          coor[3*(n10-1)+ii] = coorlocal[k][ii];
        listn[nyg*nxg*iz+nxg*iy+ix].push_back(n10);
      }
    }
  }
  n10_[0]=n10;
}
