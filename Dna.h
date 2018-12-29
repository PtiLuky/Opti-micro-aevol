//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <zlib.h>

#include "Threefry.h"

// REPETITION_MARGIN * FRAME_SIZE > 22 needed for the repetition
#define REPETITION_MARGIN 1

#define FRAME_SIZE_POW 5
#define FRAME_SIZE (1<<FRAME_SIZE_POW)
// index of the frame
#define FRAME(pos) ((pos) >> (FRAME_SIZE_POW)) // pos / FRAME_SIZE
// index in frame from left
#define INDEXFL(pos, frame) ((pos) - ((frame) << FRAME_SIZE_POW)) // pos - frame * FRAME_SIZE
// index in frame from right
#define INDEXFR(pos, frame) ((((frame) + 1) << FRAME_SIZE_POW) - ((pos) + 1)) // (frame+1) * FRAME_SIZE - (pos+1)
// left move
#define LEFT(x, a) ((x) << (a))
// right move
#define RIGHT(x, a) ((x) >> (a))
// keep a first bits on a total length of len
#define FIRST(x, a, len) ((x) >> ((len)-(a)))
// keep a last bits
#define LAST(x, a) ((x) & ((1 << (a)) - 1))
// bit at a (from right)
#define AT(x, a) (((x) >> (a)) & 0b1)
// set var count as the number of bit at 1 in x (count must be a l-value)
#define COUNT(x, count) for(count = 0 ; (x); ++count) (x) &= (x) - 1

constexpr int32_t CODON_SIZE = 3;

constexpr const uint32_t PROM_SEQ = 0b0101011001110010010110;
constexpr const int PROM_SEQ_L = 22;
constexpr const uint32_t SHINE_DAL_SEQ = 0b011011000;
constexpr const int SHINE_DAL_SEQ_L = 9;
constexpr const uint32_t PROTEIN_END = 0b001; // CODON_STOP
constexpr const int PROTEIN_END_L = 3; // CODON_STOP


class ExpManager;

struct Sequence {
  uint32_t* seq;
  int length; // size of usefull data
  int nbElem; // nb of elem
  int size; // true size
  Sequence(){
    this->length = 0;
    this->nbElem = 0;
    this->size = 0;
  }
  Sequence(int length){
    this->length = length;
    this->nbElem = ((length - 1)/(sizeof(uint32_t)*8) + 1) + REPETITION_MARGIN; // byte to bit
    this->size = this->nbElem *  sizeof(uint32_t);
    this->seq = (uint32_t*) malloc(size);
  }
  ~Sequence(){ if(length) free(seq);}
};

class Dna {

 public:
  Dna() = default;

  Dna(const Dna& clone);

  Dna(int length, Threefry::Gen& rng);

  Dna(uint32_t* genome, int length);

  Dna(int length);

  ~Dna() = default;

  int length() const;

  void save(gzFile backup_file);
  void load(gzFile backup_file);

  void set(int pos, bool c);

  void do_switch(int pos);

  int promoter_at(int pos);

  int terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

  struct Sequence seq_;
};
