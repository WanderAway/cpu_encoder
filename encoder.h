#pragma once 
#include "defines.h"
#include <vector>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <thread>

#define INTRA_H_MODE 0
#define INTRA_V_MODE 1

#define PAR_BLK 1
#define PAR_FRM 2

using namespace std;


typedef struct v2d{
  int x;
  int y; 
  v2d(int xi=0, int yi=0) {
    x = xi; y = yi;
  }
  bool operator<(v2d& rhs){if (y == rhs.y) return abs(x) < abs(rhs.x); else return abs(y) < abs(rhs.y);}
  // reverse the direction of the vector
  v2d& neg(){x = -1 * x; y = -1 * y; return *this;}
  v2d operator+(v2d& rhs){return v2d(x+rhs.x, y+rhs.y);}
  v2d operator-(v2d& rhs){return v2d(x-rhs.x, y-rhs.y);}
} vector2d;


class block;
class encoder;

class frame
{
  friend class block;
  friend class encoder;
  public: 
    frame();
    ~frame();
    frame(const frame& copy);
    void operator=(const frame& rhs);
    void swap(frame& rhs);
    
    void set_resolution(const vector2d& res);
    u8 get_px(const vector2d& px);
    void set_px(const vector2d& px, u8 val);
    u8* get_ptr(const vector2d& px);
    void set_qtc(short* qtc);
    void set_res(short* res);

    
  private: 
    unsigned size;
    vector2d resolution;
    u8* buf;
    short* qtc;
    short* res;
};

class block
{
  friend class encoder;
  public: 
    block(encoder* enc, frame* frm, unsigned short blk_size = 16, unsigned blkno = 0);
    
    void set_blk(unsigned blkno);
    u8 intra_generate(frame* targ);
    vector2d mv_generate(frame* ref, frame* targ);

    block& operator++ (int);

  private: 
    void dct_quantize(short *res);
    void dequantize_idct(short *res);
    
    
    unsigned blkno; 
    vector2d blk_loc; 
    vector2d px_loc;
    unsigned short blk_size;
    unsigned px_size;
    unsigned h_blks;
    unsigned v_blks;
    encoder* enc;
    frame* frm;
};

class encoder
{
  friend class block;
  public: 
    encoder(u8 reference_frames = 1, vector2d res = vector2d(352,288), unsigned short blk_size = 16, 
            unsigned search_range = 4, u8 parallel_mode = 0, u8 qp = 5);
    ~encoder();

    void encode(ifstream& infile, ofstream& outfile, unsigned frames = 10, unsigned i_period = 10);
    
  private: 

    void generate_iframe();
    void generate_pframe(frame* ref);

    void generate_iframe_row(unsigned row, frame* cur, frame* targ);
    void generate_pframe_row(unsigned row, frame* ref, frame* cur, frame* targ);
    
    vector<frame> ref_frms;
    frame cur_frame;
    frame cur_frame2;
    frame gen_frame;
    frame gen_frame2;
    unsigned h_blks;
    unsigned v_blks;
    u8 qp;
    u8 *i_mode;
    unsigned search_range;
    vector<vector2d> mv;
    short *residuals;
    short *residuals2;
    short *qtc;
    short *qtc2;
    short *q;
    u8 parallel_mode;
    unsigned short blk_size;
};
