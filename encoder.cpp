#include "encoder.h"

//-------------------------------------------------------------------------
// frame class 
frame::frame()
{
  buf = NULL;
}

frame::~frame()
{
  if (!buf)
    delete[] buf;
}

frame::frame(const frame& copy)
{
  this->resolution.x = copy.resolution.x;
  this->resolution.y = copy.resolution.y;
  this->size = copy.size;


  buf = new u8[size];
  memcpy(this->buf, copy.buf, copy.size);
    
}


void frame::operator=(const frame& rhs)
{
  this->resolution.x = rhs.resolution.x;
  this->resolution.y = rhs.resolution.y;
  this->size = rhs.size;

  if (this->buf)
    delete[] this->buf;

  buf = new u8[size];
  memcpy(this->buf, rhs.buf, rhs.size);
}

void frame::swap(frame& rhs)
{
  u8 *temp = this->buf;
  this->buf = rhs.buf;
  rhs.buf = temp;
}

void frame::set_resolution(const vector2d& res)
{
  // if something already here
  if (this->buf)
    delete[] this->buf;
  size = res.x * res.y;
  resolution.x = res.x;
  resolution.y = res.y;
  buf = new u8[size];
}
    

u8 frame::get_px(const vector2d& px)
{
  unsigned index = px.y * resolution.x + px.x;
  return buf[index];
}

void frame::set_px(const vector2d& px, u8 val)
{
  unsigned index = px.y * resolution.x + px.x;
  buf[index] = val;
}

u8* frame::get_ptr(const vector2d& px)
{
  unsigned index = px.y * resolution.x + px.x;
  return (buf+index);
}

void frame::set_qtc(short* qtc)
{
  this->qtc = qtc;
}

void frame::set_res(short* res)
{
  this->res = res;
}
//-------------------------------------------------------------------------
// block class 
block::block(encoder* enc, frame* frm, unsigned short blk_size, unsigned blkno)
{
  this->blkno = blkno; 
  this->blk_size = blk_size;
  this->enc = enc;
  this->px_size = blk_size * blk_size;
  this->frm = frm;
  
  h_blks = enc->h_blks;
  v_blks = enc->v_blks;

  blk_loc.y = blkno / h_blks;
  blk_loc.x = blkno % h_blks;

  px_loc.x = blk_size * blk_loc.x;
  px_loc.y = blk_size * blk_loc.y;
  
  // ideally this wouldn't happen. 
  if (blk_loc.y > v_blks) 
    cerr << "WARNING: block " << blkno << " out of range\n";
}

void block::set_blk(unsigned blkno)
{
  this->blkno = blkno;
  blk_loc.y = blkno / h_blks;
  blk_loc.x = blkno % h_blks;
  px_loc.x = blk_size * blk_loc.x;
  px_loc.y = blk_size * blk_loc.y;
}

u8 block::intra_generate(frame* targ)
{
  vector2d rpx;
  u8 mode = INTRA_H_MODE;
  u8 *hpb, *vpb;
  hpb = new u8[px_size];
  vpb = new u8[px_size];

  
  u8 hlen;
  // find number of px left in block/frame horizontally
  if (blk_size + px_loc.x > frm->resolution.x)
    hlen = frm->resolution.x - px_loc.x;
  else 
    hlen = blk_size;

  u8 vlen; 
  // same but vertically
  if (blk_size + px_loc.y > frm->resolution.y)
    vlen = frm->resolution.y - px_loc.y;
  else 
    vlen = blk_size;

  // do horizontal prediction:
  if (px_loc.x == 0) // leftmost block
    memset(hpb, 128, px_size);
  else 
  {
    rpx.x = px_loc.x - 1;
    for (unsigned i=0; i<vlen; i++)
    {
      rpx.y = i + px_loc.y;
      u8 temp = targ->get_px(rpx);
      memset(hpb, temp, hlen);
    }
  }

  // vertical prediction
  if (px_loc.y == 0) // topmost block
    memset(vpb, 128, px_size);
  else 
  {
    rpx.y = px_loc.y - 1;
    for (unsigned i=0; i<hlen; i++)
    {
      rpx.x = i + px_loc.x;
      u8 temp = targ->get_px(rpx);
      unsigned index = i;

      for (unsigned j=0; j<vlen; j++)
      {
        vpb[index] = temp;
        index += blk_size;
      }
    }
  }

  
  // calculate SAD for both
  short hres[px_size];
  short vres[px_size];

  unsigned hsad=0, vsad=0;
  unsigned pxno;
  for (unsigned i=0; i<vlen; i++)
  {
    rpx.y = px_loc.y + i;
    pxno = i * blk_size;
    for (unsigned j=0; j<hlen; j++)
    {
      rpx.x = px_loc.x + j;
      short temp = frm->get_px(rpx);
      hres[pxno] = temp - hpb[pxno];
      vres[pxno] = temp - vpb[pxno];
      hsad += abs(hres[pxno]);
      vsad += abs(vres[pxno]);

      pxno++;
    }
  }

  u8 *best;
  // save one with least SAD
  if (vsad < hsad)
  {
    enc->i_mode[blkno] = INTRA_V_MODE;
    dct_quantize(vres);
    best = vpb;
    mode = INTRA_V_MODE;
  }
  else 
  {
    enc->i_mode[blkno] = INTRA_H_MODE;
    dct_quantize(hres);
    best = hpb;
    mode = INTRA_H_MODE;
  }
  short *restarg = frm->res + blkno * px_size;
  dequantize_idct(restarg);

  unsigned index;
  // write block to generated buffer. 
  for (unsigned i=0; i<vlen; i++)
  {
    index = i * blk_size;
    rpx.y = px_loc.y + i;
    for (unsigned j=0; j<hlen; j++)
    {
      rpx.x = px_loc.x + j;
      short temp;
      temp = best[index] + restarg[index];
      if (temp < 0)
        temp = 0;
      else if (temp > 255)
        temp = 255;
      
      targ->set_px(rpx, temp);

      index++;
    }
  }

  enc->i_mode[blkno] = mode;
  delete[] vpb;
  delete[] hpb;
  return mode;
}


vector2d block::mv_generate(frame* ref, frame* targ)
{
  vector2d mv;
  int r = enc->search_range;
  int lbound, rbound, tbound, bbound;
  
  // define search boundry
  if (px_loc.x - r < 0)
    lbound = 0;
  else 
    lbound = px_loc.x - r;
  
  if (px_loc.x + r + blk_size >= frm->resolution.x)
    rbound = frm->resolution.x - blk_size;
  else 
    rbound = px_loc.x + r;
  
  if (px_loc.y - r < 0)
    tbound = 0;
  else 
    tbound = px_loc.y - 4;

  if (px_loc.y + blk_size + r > frm->resolution.y)
    bbound = frm->resolution.y - blk_size;
  else 
    bbound = px_loc.y + r;
  
  // handle edge case: 
  if (lbound > rbound || tbound > bbound) // block out of bound
  { // TODO not within our use case yet
    
  }
  else 
  {
    vector2d ofst;
    vector2d rpx, cpx;
    unsigned sad, lsad=0; 
    bool first = true;
    unsigned index;
    short res[px_size], lres[px_size];


    u8 hlen;
    // find number of px left in block/frame horizontally
    if (blk_size + px_loc.x > frm->resolution.x)
      hlen = frm->resolution.x - px_loc.x;
    else 
      hlen = blk_size;

    u8 vlen; 
    // same but vertically
    if (blk_size + px_loc.y > frm->resolution.y)
      vlen = frm->resolution.y - px_loc.y;
    else 
      vlen = blk_size;
    for (unsigned i=tbound; i<=bbound; i++)
    {
      ofst.y = i - px_loc.y;
      for (unsigned j=lbound; j<=rbound; j++)
      {
        ofst.x = j - px_loc.x;
        // calculate sad here
        cpx = px_loc;
        sad=0;
        index=0;
        for (rpx.y = i; rpx.y<i+blk_size/* && rpx.y<ref->resolution.y*/; rpx.y++)
        {
          for (rpx.x = j; rpx.x<j+blk_size/* && rpx.x<ref->resolution.x*/; rpx.x++)
          {
            res[index] = (short)(frm->get_px(cpx) ) - ref->get_px(rpx);
            sad += abs(res[index]);
            index++;
            cpx.x++;
          }
          cpx.x = px_loc.x;
          cpx.y++;
        }

        // save best
        if (first || sad < lsad)
        {
          lsad = sad;
          mv = ofst;
          memcpy(lres, res, sizeof(short)*px_size);
          first = false;
        }
      }
    }
    // store best in residuals after dct/quantize, generate frame. 
    dct_quantize(lres);
    short *restarg = frm->res + blkno * px_size;
    dequantize_idct(restarg);

    
    // write block to generated buffer. 
    rpx = mv + px_loc;
    cpx = px_loc;
    int rx = rpx.x;
    int cx = cpx.x;
    index = 0;
    for (unsigned i=0; i<vlen; i++)
    {
      for (unsigned j=0; j<hlen; j++)
      {
        short temp;
        temp = ref->get_px(rpx) + restarg[index]; 
        if (temp < 0)
          temp = 0;
        else if (temp > 255)
          temp = 255;
        
        targ->set_px(cpx, temp);

        rpx.x++;
        cpx.x++;
        index++;
      }
      rpx.x = rx;
      cpx.x = cx;
      rpx.y++;
      cpx.y++;
    }

  }
  
  enc->mv[blkno] = mv;
  return mv;
}


block& block::operator++(int)
{
  blkno++;
  blk_loc.x += 1;
  if (blk_loc.x >= h_blks)
  {
    blk_loc.x = 0;
    blk_loc.y += 1;
  }

  px_loc.x = blk_size * blk_loc.x;
  px_loc.y = blk_size * blk_loc.y;

  if (blk_loc.y > v_blks) 
    cerr << "WARNING: block " << blkno << " out of range\n";
  return *this;
}


void block::dct_quantize(short *res)
{
  float ci, cj, dct1, sum;
  short coef[px_size];
  
  for (unsigned i=0; i<blk_size; i++)
  {
    if (i==0)
      ci = 1 / sqrt(blk_size);
    else ci = sqrt(2) / sqrt(blk_size);
    
    for (unsigned j=0; j<blk_size; j++)
    {
      if (j == 0)
        cj = 1/sqrt(blk_size);
      else 
        cj = sqrt(2) / sqrt(blk_size);

      sum = 0;
      
      for (unsigned k=0; k<blk_size; k++)
      {
        for (unsigned l=0; l < blk_size; l++)
        {
          dct1 = res[k * blk_size + l] * 
                 cos(((2*k + 1) * i * PI) / (2 * blk_size)) * 
                 cos(((2*l + 1) * j * PI) / (2 * blk_size));
          sum += dct1;
        }
      }
      coef[i*blk_size + j] = round(ci * cj * sum);
    }
  }
  
  short *bq = frm->qtc + blkno * px_size;
  unsigned bpno=0;

  for (unsigned i=0; i<blk_size; i++)
  {
    for (unsigned j=0; j<blk_size; j++)
    {
      short temp; 
      temp = round( ( (float)coef[bpno] / enc->q[bpno]));
      bq[bpno] = temp;
      bpno++;
    }
  }
}

void block::dequantize_idct(short *res)
{
  short *bq = frm->qtc + blkno * px_size;
  short coef[px_size];
  unsigned bpno = 0;

  for (unsigned i=0; i<blk_size; i++)
  {
    for (unsigned j=0; j<blk_size; j++)
    {
      short temp; 
      temp = round(bq[bpno] * enc->q[bpno]);
      coef[bpno] = temp;
      bpno++;
    }
  }

  float intermediate[px_size];
  
  // horizontal idct
  for (unsigned i=0; i<blk_size; i++)
  {
    for(unsigned j=0; j<blk_size; j++)
    {
      float val=0.0;
      for (unsigned k=0; k<blk_size; k++)
      {
        if (k==0)
          val += coef[i*blk_size] / sqrt(2);
        else 
          val += coef[i*blk_size+k] * cos( PI * k * (j+0.5) / blk_size);
      }
      val *= sqrt(2.0/blk_size);
      intermediate[i*blk_size + j] = val;
    }
  }

  for (unsigned i=0; i<blk_size; i++)
  {
    for (unsigned j=0; j<blk_size; j++)
    {
      
      float val=0.0;
      for (unsigned k=0; k<blk_size; k++)
      {
        if (k == 0)
          val += intermediate[i] / sqrt(2);
        else 
          val += intermediate[k*blk_size+i] * cos( PI * k * (j+0.5) / blk_size);
      }
      val *= sqrt(2.0/blk_size);
      res[j*blk_size+i] = round(val);
    }
  }
}


//-------------------------------------------------------------------------
// encoder class 
encoder::encoder(u8 reference_frames, vector2d res, unsigned short blk_size, unsigned search_range, u8 parallel_mode, u8 qp)
{
  ref_frms.resize(reference_frames);
  
  ref_frms[0].set_resolution(res);

  cur_frame.set_resolution(res);
  cur_frame2.set_resolution(res);
  
  gen_frame.set_resolution(res);
  gen_frame2.set_resolution(res);

  this->blk_size = blk_size;
  this->residuals = new short[cur_frame.size];
  this->residuals2 = new short[cur_frame.size];
  this->qtc = new short[cur_frame.size];
  this->qtc2 = new short[cur_frame.size];
  this->parallel_mode = parallel_mode;
  this->qp = qp;
  this->search_range = search_range;

  cur_frame.set_qtc(qtc);
  cur_frame.set_res(residuals);
  cur_frame2.set_res(residuals2);
  cur_frame2.set_qtc(qtc2);

  this->q = new short[blk_size*blk_size];

  for (unsigned i=0; i<blk_size; i++)
  {
    for (unsigned j=0; j<blk_size; j++)
    {
      if ( (i+j) < blk_size - 1)
        q[i*blk_size + j] = 1 << qp;
      else if ( (i+j) == blk_size)
        q[i*blk_size+j] = 1 << (qp+1);
      else 
        q[i*blk_size+j] = 1 << (qp+2);
    }
  }
  
  h_blks = cur_frame.resolution.x / blk_size;
  if (cur_frame.resolution.x % blk_size)
    h_blks++;
  v_blks = cur_frame.resolution.y / blk_size;
  if (cur_frame.resolution.y % blk_size)
    v_blks++;

  unsigned blk_cnt = h_blks * v_blks;

  i_mode = new u8[blk_cnt];
  mv.resize(h_blks*v_blks);

}

encoder::~encoder()
{
  delete[] residuals;
  delete[] qtc;
  delete[] qtc2;
  delete[] q;
  delete[] i_mode;
}


void encoder::encode(ifstream& infile, ofstream& outfile, unsigned frames, unsigned i_period)
{
  unsigned count=0;
  if (parallel_mode == PAR_FRM)
  {
    thread par_thread;
    while (count < frames)
    {
      cout << "frame " << count << endl;
      cout << "frame " << count+1 << endl;
      infile.read((char*)cur_frame.buf, cur_frame.size);
      if (count != frames-1)
        infile.read((char*)cur_frame2.buf, cur_frame2.size);
      if (count % i_period)
      {
        // generate first two rows
        generate_pframe_row(0, &ref_frms[0], &cur_frame, &gen_frame);
        generate_pframe_row(1, &ref_frms[0], &cur_frame, &gen_frame);
      }
      else 
      {
        generate_iframe_row(0, &cur_frame, &gen_frame);
        generate_iframe_row(1, &cur_frame, &gen_frame);
      }
      unsigned row2 = 0;
      
      for (unsigned i=2; i<v_blks; i++)
      {
        if (par_thread.joinable())
        {
          par_thread.join();
        }

        // start thread: 
        if (count != frames-1)
        {
          par_thread = thread(&encoder::generate_pframe_row, this, row2, &gen_frame, &cur_frame2, &gen_frame2);
          row2++;
        }

        // generate row on this thread 
        if (count % i_period)
          generate_pframe_row(i, &ref_frms[0], &cur_frame, &gen_frame);
        else 
          generate_iframe_row(i, &cur_frame, &gen_frame);
      }
      if (par_thread.joinable())
      {
        par_thread.join();
      }

      // finish processing the second frame
      while (row2 < v_blks && count != frames-1)
      {
        generate_pframe_row(row2, &gen_frame, &cur_frame2, &gen_frame2);
        row2++;
      }


      // write to stream: 
      outfile.write((char*)gen_frame.buf, gen_frame.size);
      outfile.write((char*)gen_frame2.buf, gen_frame2.size);

      // use gen_frame 2 as reference
      ref_frms[0].swap(gen_frame2);

      count += 2;
    }

  }
  else 
  {
    frame* ref = &ref_frms[0];
    while (count < frames) 
    {
      cout << "frame " << count << endl;
      // read from infile to cur_frame
      infile.read((char*)cur_frame.buf, cur_frame.size);

      if (count % i_period)
      {// p frame
        generate_pframe(ref);
      }
      else 
      {// i frame
        generate_iframe();
      }

      // write gen frame to output 
      outfile.write((char*)gen_frame.buf, gen_frame.size);

      // swap gen frame and ref frame so that generated becomes previous
      ref->swap(gen_frame);
      count++;
    }
  }
}

// block level parallism always enabled now
void encoder::generate_iframe()
{
  // initialize second to first block on second row
  block first(this, &cur_frame), second(this, &cur_frame, 16, h_blks);
  unsigned blkno=0;
  unsigned row = 0;
  bool gen_thread = true;
  thread blk_thrd;

  if (parallel_mode == 1)
  { // block level parallization
    while (row < v_blks)
    {
      // calculate blkno
      blkno = row * h_blks;
      first.set_blk(blkno);
      // don't generate next row when already last row
      if (row == v_blks-1)
      {
        gen_thread = false;
        if (blk_thrd.joinable())
          blk_thrd.join(); 
      }
      else
      {
        // update second after sync
        if (blk_thrd.joinable())
          blk_thrd.join(); 
        second.set_blk(blkno + h_blks);
      }
      
      // generate all blocks in row
      for (unsigned i=0; i<h_blks; i++)
      {
        first.intra_generate(&gen_frame); 

        if (gen_thread)
        {
          // synchronize
          if (blk_thrd.joinable())
          {
            blk_thrd.join(); 
            second++;
          }
          // once the block is done, the one below it can go ahead. 
          blk_thrd = thread(&block::intra_generate, &second, &gen_frame);
        }
        first++;
      }

      row += 2;
    }
    // sync one last time before exiting
    if (blk_thrd.joinable())
      blk_thrd.join();
  }
  else 
  {
    for (unsigned i=0; i<h_blks * v_blks; i++)
    {
      first.intra_generate(&gen_frame);
      first++;
    }
  }

}


void encoder::generate_pframe(frame* ref)
{
  // initialize second to first block on second row
  block first(this, &cur_frame), second(this, &cur_frame, 16, h_blks);
  unsigned blkno=0;
  unsigned row = 0;
  bool gen_thread = true;
  thread blk_thrd;

  if (parallel_mode == 1 || parallel_mode == 3)
  { // block level parallization
    while (row < v_blks)
    {
      // calculate blkno
      blkno = row * h_blks;
      first.set_blk(blkno);
      // don't generate next row when already last row
      if (row == v_blks-1)
      {
        gen_thread = false;
        if (blk_thrd.joinable())
          blk_thrd.join(); 
      }
      else
      {
        // update second after sync
        if (blk_thrd.joinable())
          blk_thrd.join(); 
        second.set_blk(blkno + h_blks);
      }
      
      // generate all blocks in row
      for (unsigned i=0; i<h_blks; i++)
      {
        first.mv_generate(ref, &gen_frame); 

        if (gen_thread)
        {
          // synchronize
          if (blk_thrd.joinable())
          {
            blk_thrd.join(); 
            second++;
          }
          // once the block is done, the one below it can go ahead. 
          blk_thrd = thread(&block::mv_generate, &second, ref, &gen_frame);
        }
        first++;
      }

      row += 2;
    }
    // sync one last time before exiting
    if (blk_thrd.joinable())
      blk_thrd.join();
  }
  else 
  {
    for (unsigned i=0; i<h_blks * v_blks; i++)
    {
      first.mv_generate(ref, &gen_frame);
      first++;
    }
  }
  
}


void encoder::generate_iframe_row(unsigned row, frame* cur, frame* targ)
{
  unsigned blkno = row * h_blks;
  block blk(this, cur, 16, blkno);
  
  for (unsigned i=0; i<h_blks; i++)
  {
    blk.intra_generate(targ);
    blk++;
  }

}

void encoder::generate_pframe_row(unsigned row, frame* ref, frame* cur, frame* targ)
{
  unsigned blkno = row * h_blks;
  block blk(this, cur, 16, blkno);
  
  for (unsigned i=0; i<h_blks; i++)
  {
    blk.mv_generate(ref, targ);
    blk++;
  }
}
