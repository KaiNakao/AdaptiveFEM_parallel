#include "hdf5.h"
#include "hdf5_hl.h"
#include "stdlib.h"
#include <omp.h>

#define I32 H5T_NATIVE_INT
#define I64 H5T_NATIVE_LONG
#define R64 H5T_NATIVE_DOUBLE

//#define PRTDBG
//#define PRTPERF
double prfthresh = 5.0;

#define GLBHDFSIZE 100
hid_t glb_file_id[GLBHDFSIZE];
int glb_file_id_open[GLBHDFSIZE];
int init = 0;

static herr_t
my_H5LT_make_dataset_numerical( hid_t loc_id, const char *dset_name, int rank, const hsize_t *dims, hid_t tid,  const void *data )
{
  hid_t   did = -1, sid = -1;

  /* Create the data space for the dataset. */
  if((sid = H5Screate_simple(rank, dims, NULL)) < 0)
    return -1;

  /* Create the dataset. */
  if((did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
    goto out;

  /* Write the dataset only if there is data to write */
  if(data)
    if(H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
      goto out;

  /* End access to the dataset and release resources used by it. */
  if(H5Dclose(did) < 0)
    return -1;

  /* Terminate access to the data space. */
  if(H5Sclose(sid) < 0)
    return -1;
  
  return 0;

out:
  H5E_BEGIN_TRY {
    H5Dclose(did);
    H5Sclose(sid);
  } H5E_END_TRY;
  return -1;
}

static herr_t
my_H5LT_read_dataset_numerical(hid_t loc_id, const char *dset_name, hid_t tid, void *data)
{
  hid_t   did;

  /* Open the dataset. */
  if((did = H5Dopen2(loc_id, dset_name, H5P_DEFAULT)) < 0)
    return -1;

  /* Read */
  if(H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
    goto out;

  /* End access to the dataset and release resources used by it. */
  if(H5Dclose(did))
    return -1;

  return 0;

out:
  H5Dclose(did);
  return -1;
}


void set_null_char_( char* filename, long* pos )
{
  filename[*pos] = '\0';
}

/* initialize file ids */
void hdf_initialize()
{
  int i;
  for( i = 0; i < GLBHDFSIZE; ++i ) glb_file_id_open[i] = 0;
}

/* create or open file */
void hdf_create_open_file( const char* filename, long* fid, int type )
{
  if( init == 0 )
  {
    init = 1;
    hdf_initialize();
  }
  if( *fid < 0 || *fid >= GLBHDFSIZE )
  {
    printf("open_hdf_file: fid size is larger than GLBHDFSIZE\n");
    fflush(stdout);
    exit(-1);
  }
  if( glb_file_id_open[*fid] != 0 )
  {
    printf("open_hdf_file: fid is already used\n");
    fflush(stdout);
    exit(-1);
  }
  glb_file_id_open[*fid] = 1;
  if(type==0)
    glb_file_id[*fid] = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if(type==1)
    glb_file_id[*fid] = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
}

void hdf_create_file_( const char* filename, long* fid )
{
  hdf_create_open_file( filename, fid, 0 );
}

void hdf_open_file_( const char* filename, long* fid )
{
  hdf_create_open_file( filename, fid, 1 );
}

/* write double array */
void hdf_write_double_array_( long* fid, char* type, double* array, long* size )
{
  double t;
  int rank;
  hsize_t dims[1];
  herr_t status;

  if( glb_file_id_open[*fid] == 0 )
  {
    printf("file is not open\n");
    fflush(stdout);
    exit(-1);
  }
  
  t = omp_get_wtime();
  dims[0]=*size;
  rank = 1;
  status = my_H5LT_make_dataset_numerical(glb_file_id[*fid],type,rank,dims,R64,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" write %s took (s): %f\n",type,t);
  fflush(stdout);
#endif
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" write info of double array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*8);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*8.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* write long array */
void hdf_write_long_array_( long* fid, char* type, long* array, long* size )
{
  double t;
  int rank;
  hsize_t dims[1];
  herr_t status;

  if( glb_file_id_open[*fid] == 0 )
  {
    printf("file is not open\n");
    fflush(stdout);
    exit(-1);
  }
  
  t = omp_get_wtime();
  dims[0]=*size;
  rank = 1;
  status = my_H5LT_make_dataset_numerical(glb_file_id[*fid],type,rank,dims,I64,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" write %s took (s): %f\n",type,t);
  fflush(stdout);
#endif 
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" write info of long array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*8);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*8.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* write int array */
void hdf_write_int_array_( long* fid, char* type, int* array, long* size )
{
  double t;
  int rank;
  hsize_t dims[1];
  herr_t status;

  if( glb_file_id_open[*fid] == 0 )
  {
    printf("file is not open\n");
    fflush(stdout);
    exit(-1);
  }
  
  t = omp_get_wtime();
  dims[0]=*size;
  rank = 1;
  status = my_H5LT_make_dataset_numerical(glb_file_id[*fid],type,rank,dims,I32,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" write %s took (s): %f\n",type,t);
  fflush(stdout);
#endif 
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" write info of long array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*4);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*4.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* check dataset existence (1: exists, 0: not exists) */
int hdf_has_dataset_( long* fid, char* type )
{
  if( *fid < 0 || *fid >= GLBHDFSIZE ) return 0;
  if( !glb_file_id_open[*fid] ) return 0;
  htri_t exists = H5Lexists(glb_file_id[*fid], type, H5P_DEFAULT);
  return (exists > 0) ? 1 : 0;
}

/* get array size */
void hdf_get_array_size_( long* fid, char* type, long* size )
{
  hsize_t dims[1];
  herr_t status;

  status = H5LTget_dataset_info(glb_file_id[*fid],type,dims,NULL,NULL);
  *size = dims[0];
}

/* read double array */
void hdf_read_double_array_( long* fid, char* type, double* array, long* memsize, long* size )
{
  double t;
  hsize_t dims[1];
  herr_t status;

  t = omp_get_wtime();
  status = H5LTget_dataset_info(glb_file_id[*fid],type,dims,NULL,NULL);
  *size = dims[0];
  if(*memsize < *size )
  {
    printf("not enough memory in %s %ld %ld",type,*memsize,*size);
    fflush(stdout);
    exit(0);
  }
  status = my_H5LT_read_dataset_numerical(glb_file_id[*fid],type,R64,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" read %s took (s): %f\n",type,t); 
  fflush(stdout);
#endif 
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" read info of double array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*8);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*8.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* read long array */
void hdf_read_long_array_( long* fid, char* type, long* array, long* memsize, long* size )
{
  double t;
  hsize_t dims[1];
  herr_t status;
  long i;

  t = omp_get_wtime();
  status = H5LTget_dataset_info(glb_file_id[*fid],type,dims,NULL,NULL);
  *size = dims[0];
  if(*memsize < *size )
  {
    printf("not enough memory in %s %ld %ld",type,*memsize,*size);
    fflush(stdout);
    exit(0);
  }
  status = my_H5LT_read_dataset_numerical(glb_file_id[*fid],type,I64,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" read %s took (s): %f\n",type,t);
  printf(" status ld %ld\n",status);
  printf(" status d %d\n",status);
  for( i = 0; i < 10; ++i )
  {
    if( i == *size )
      break;
    printf("val[%ld]: %ld\n",i,array[i]);
  }
  fflush(stdout);
#endif 
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" read info of long array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*8);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*8.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* read int array */
void hdf_read_int_array_( long* fid, char* type, int* array, long* memsize, long* size )
{
  double t;
  hsize_t dims[1];
  herr_t status;
  long i;

  t = omp_get_wtime();
  status = H5LTget_dataset_info(glb_file_id[*fid],type,dims,NULL,NULL);
  *size = dims[0];
  if(*memsize < *size )
  {
    printf("not enough memory in %s %ld %ld",type,*memsize,*size);
    fflush(stdout);
    exit(0);
  }
  status = my_H5LT_read_dataset_numerical(glb_file_id[*fid],type,I32,array);
  t = omp_get_wtime() - t;
#ifdef PRTDBG
  printf(" read %s took (s): %f\n",type,t);
  printf(" status ld %ld\n",status);
  printf(" status d %d\n",status);
  for( i = 0; i < 10; ++i )
  {
    if( i == *size )
      break;
    printf("val[%ld]: %ld\n",i,array[i]);
  }
  fflush(stdout);
#endif 
#ifdef PRTPERF
  if( t > prfthresh )
  {
    printf("**************************\n");
    printf(" read info of long array: %s\n",type);
    printf(" size (B): %ld\n",(*size)*4);
    printf(" time (s): %f\n",t);
    printf(" rate (B/s): %f\n",(*size)*4.0/t);
    printf("**************************\n");
    fflush(stdout);
  }
#endif
}

/* get long array component */
void hdf_read_long_array_part_( long* fid, char* type, long* start, long* size, long* val )
{
   hsize_t dims[1], dimsm[1];
   hid_t dataset_id;
   hid_t       dataspace_id, memspace_id;
   herr_t      status;
   hsize_t count[1];              /* size of subset in the file */
   hsize_t offset[1];             /* subset offset in the file */
   hsize_t stride[1];
   hsize_t block[1];
   long i, j;
   int rank;

   rank = 1;
   dataset_id = H5Dopen2 (glb_file_id[*fid], type, H5P_DEFAULT);
   offset[0] = *start;
   count[0] = *size;
   stride[0] = 1;
   block[0] = 1;

   dimsm[0] = *size;
   memspace_id = H5Screate_simple (rank, dimsm, NULL);
   dataspace_id = H5Dget_space (dataset_id);
   status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, stride, count, block);

   status = H5Dread (dataset_id, I64, memspace_id, dataspace_id, H5P_DEFAULT, val);
   status = H5Sclose (memspace_id);
   status = H5Sclose (dataspace_id);
   status = H5Dclose (dataset_id);
}

/* get int array component */
void hdf_read_int_array_part_( long* fid, char* type, long* start, long* size, int* val )
{
   hsize_t dims[1], dimsm[1];
   hid_t dataset_id;
   hid_t       dataspace_id, memspace_id;
   herr_t      status;
   hsize_t count[1];              /* size of subset in the file */
   hsize_t offset[1];             /* subset offset in the file */
   hsize_t stride[1];
   hsize_t block[1];
   long i, j;
   int rank;

   rank = 1;
   dataset_id = H5Dopen2 (glb_file_id[*fid], type, H5P_DEFAULT);
   offset[0] = *start;
   count[0] = *size;
   stride[0] = 1;
   block[0] = 1;

   dimsm[0] = *size;
   memspace_id = H5Screate_simple (rank, dimsm, NULL);
   dataspace_id = H5Dget_space (dataset_id);
   status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, stride, count, block);

   status = H5Dread (dataset_id, I32, memspace_id, dataspace_id, H5P_DEFAULT, val);
   status = H5Sclose (memspace_id);
   status = H5Sclose (dataspace_id);
   status = H5Dclose (dataset_id);
}

/* close file */
void hdf_close_file_( long* fid )
{
  herr_t status;
  if( glb_file_id_open[*fid] == 0 )
  {
    printf("close_hdf_file: fid is already closed\n");
    fflush(stdout);
    exit(-1);
  }
  status = H5Fclose(glb_file_id[*fid]);
  glb_file_id_open[*fid] = 0;
}
