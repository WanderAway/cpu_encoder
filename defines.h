#pragma once
#include <cstdint>

#define u8 uint8_t 
#define u16 uint16_t
#define u32 uint32_t
#define i8 int8_t
#define i16 int16_t
#define i32 int32_t

typedef struct stat STAT;
#define PI 3.14159265358979323846

//---------------------------------------------------------
// output options
//#define MV_TOFILE 1
//#define MODE_TOFILE 1
//#define RES_TOFILE 1
//#define ENTROPY_TOFILE 1

#define DEBUG
//---------------------------------------------------------


#ifndef MV_TOFILE
#define MV_TOFILE 0
#endif
#ifndef MODE_TOFILE
#define MODE_TOFILE 0
#endif
#ifndef RES_TOFILE
#define RES_TOFILE 0
#endif
#ifndef ENTROPY_TOFILE
#define ENTROPY_TOFILE 0
#endif

