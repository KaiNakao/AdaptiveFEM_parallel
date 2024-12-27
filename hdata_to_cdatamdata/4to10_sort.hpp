#include <list>

bool near( std::list<coor_pair>::iterator a, std::list<coor_pair>::iterator b)
{
  double x(a->x[0]-b->x[0]);
  double y(a->x[1]-b->x[1]);
  double z(a->x[2]-b->x[2]);
  return (x*x+y*y+z*z) < global_coor_pair_eps*global_coor_pair_eps;
}

void convert_4to10_sort_cpp( long* n4_, long* ne_, long* coorsize_, double* coor, long* conn,
long* nl_, double* dsg_,long* nxg_, long* nyg_,long* nzg_,double* co,long* n10_,long* flag_)
{
  global_coor_pair_eps = 0.001;
  long n4,ne,nxg,nyg,nzg,n10,coorsize;
  double dsg;
  long i,j,ii,k;
  long cnyele[10][2];
  double coorlocal[10][3];
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
  
  std::list<coor_pair> pairs;
  
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
    
    for( k = 4; k < 10; ++k )
      pairs.push_back(coor_pair(11*i+k,coorlocal[k][0],coorlocal[k][1],coorlocal[k][2],-1));
  }
  
  pairs.sort();
  
  std::list<coor_pair>::iterator it;
  std::list<coor_pair>::iterator it1;
 
  for( it = pairs.begin(); it != pairs.end(); ++it1 )
  {
    n10++;
    for( it1 = it; it1 != pairs.end() && near(it,it1); ++it1 )
    {
      it1->type = n10;
    }
    if( (n10+2) == coorsize )
    {
      std::cout << "coorsize too small " << coorsize << "\n";
      exit(0);
    }
    if(it1 == pairs.end()) break;
    it=it1;
  }
  
  for( it = pairs.begin(); it != pairs.end(); ++it )
  {
    if(it->id == -1)
      Error("something wrong in 4to10 sort pairs");
    conn[it->id] = it->type;
  }
  
  
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
    for( k = 4; k < 10; ++k )
      for( ii = 0; ii < 3; ++ii )
        coor[3*(conn[11*i+k]-1)+ii] = coorlocal[k][ii];
  }
  
  n10_[0]=n10;
}
