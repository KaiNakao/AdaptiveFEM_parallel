#ifndef ANALYSIS_INTEGRATED_COMMON_HDF_LIB_H_
#define ANALYSIS_INTEGRATED_COMMON_HDF_LIB_H_

#ifdef __cplusplus
extern "C" {
#endif

void set_null_char_(char* filename, long* pos);

void hdf_initialize(void);
void hdf_create_open_file(const char* filename, long* fid, int type);
void hdf_create_file_(const char* filename, long* fid);
void hdf_open_file_(const char* filename, long* fid);

void hdf_write_double_array_(long* fid, char* type, double* array, long* size);
void hdf_write_long_array_(long* fid, char* type, long* array, long* size);
void hdf_write_int_array_(long* fid, char* type, int* array, long* size);

void hdf_get_array_size_(long* fid, char* type, long* size);

void hdf_read_double_array_(long* fid, char* type, double* array, long* memsize, long* size);
void hdf_read_long_array_(long* fid, char* type, long* array, long* memsize, long* size);
void hdf_read_int_array_(long* fid, char* type, int* array, long* memsize, long* size);

void hdf_read_long_array_part_(long* fid, char* type, long* start, long* size, long* val);
void hdf_read_int_array_part_(long* fid, char* type, long* start, long* size, int* val);

void hdf_close_file_(long* fid);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // ANALYSIS_INTEGRATED_COMMON_HDF_LIB_H_
