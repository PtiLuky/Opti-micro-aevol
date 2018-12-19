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

#define FRAME_SIZE_POW 3
#define FRAME_SIZE (1<<3)
// index of the frame
#define FRAME(pos) ((pos) >> (FRAME_SIZE_POW)) // pos / 8
// index in frame from left
#define INDEXFL(pos, frame) ((pos) - ((frame) << FRAME_SIZE_POW)); // pos - frame * 8
// index in frame from right
#define INDEXFR(pos, frame) ((((frame) + 1) << FRAME_SIZE_POW) - ((pos) + 1)); // (frame+1) * 8 - (pos+1)

constexpr int8_t CODON_SIZE = 3;

constexpr const uint32_t PROM_SEQ = 0b0101011001110010010110;
constexpr const int PROM_SEQ_L = 22;
constexpr const uint32_t SHINE_DAL_SEQ = 0b011011000;
constexpr const int SHINE_DAL_SEQ_L = 9;
constexpr const uint32_t PROTEIN_END = 0b001; // CODON_STOP
constexpr const int PROTEIN_END_L = 3; // CODON_STOP


class ExpManager;

struct Sequence {
  uint8_t* seq;
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
    this->nbElem = ((length - 1)/sizeof(uint8_t) + 1);
    this->size = this->nbElem *  sizeof(uint8_t);
    this->seq = (uint8_t*) malloc(size);
  }
};

class Dna {

 public:
  Dna() = default;

  Dna(const Dna& clone);

  Dna(int length, Threefry::Gen& rng);

  Dna(uint8_t* genome, int length);

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

  // std::vector<char> seq_;
  struct Sequence seq_;
};
