//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_.length){ 
  memcpy(seq_.seq, clone.seq_.seq, seq_.size);
}

Dna::Dna(int length, Threefry::Gen& rng) : seq_(length) {
  int8_t temp[length];
  
  // Generate a random genome
  for (int32_t i = 0; i < length; i++) {
    temp[i] = rng.random(NB_BASE);
  }
  for (int32_t i = 0; i < seq_.nbElem - REPETITION_MARGIN; i++) {
    seq_.seq[i] = 0;
    for(int32_t j = 0; j < FRAME_SIZE; ++j){
      seq_.seq[i] += temp[i * FRAME_SIZE + j] << (FRAME_SIZE - j - 1);
    }
  }

  /*
  // faster generation, but it changes the seed : 
  // Generate a random genome
  int32_t index_bit = 0;
  for (int32_t i = 0; i < seq_.nbElem - REPETITION_MARGIN; i++) {
    seq_.seq[i] = rng.random(UINT8_MAX);
  }
  
*/
  // repeat
  int offset = length - ((length >> FRAME_SIZE_POW) << FRAME_SIZE_POW);
  // TODO adapter la répétition pour les non-multiples de 8 !
  if(offset != 0)
    std::cout << "Please, the genome size must be a multiple of 8 (adaptation to be done)" << std::endl;
  
  for (int32_t i = seq_.nbElem - REPETITION_MARGIN; i < seq_.nbElem; i++ ){
    seq_.seq[i] = seq_.seq[i - (seq_.nbElem - REPETITION_MARGIN)];
  }
}

Dna::Dna(uint8_t* genome, int length) : seq_(length) {
  memcpy(seq_.seq, genome, seq_.size);
}

Dna::Dna(int length) : seq_(length) {
}

int Dna::length() const {
  return seq_.length;
}

void Dna::save(gzFile backup_file) {
  int dna_length = length();
  int dna_size = seq_.size;
  gzwrite(backup_file, &dna_length, sizeof(dna_length));
  gzwrite(backup_file, &dna_size, sizeof(dna_size));
  gzwrite(backup_file, seq_.seq, seq_.size);
}

void Dna::load(gzFile backup_file) {
  int dna_length;
  int dna_size;
  gzread(backup_file,&dna_length,sizeof(dna_length));
  gzread(backup_file,&dna_size,sizeof(dna_size));
  gzread(backup_file, seq_.seq, seq_.size);
}

void Dna::set(int pos, bool c) {
  // those macro are defined in Dna.h
  int frame = FRAME(pos);
  int p_in_frame_from_r = INDEXFR(pos, frame);
  uint8_t mask = 1 << p_in_frame_from_r;

  if(c) seq_.seq[frame] |= mask; // set 1
  else seq_.seq[frame] &= ~mask; // set 0
}

void Dna::do_switch(int pos) {
  int frame = FRAME(pos);
  int p_in_frame_from_r = INDEXFR(pos, frame);
  uint8_t mask = 1 << p_in_frame_from_r;

  seq_.seq[frame] ^= mask;
}


int Dna::promoter_at(int pos) {
  // TODO modulo et sortie de range...
  int frame = FRAME(pos);
  int nb_elem_in_first_frame = FRAME_SIZE - INDEXFL(pos, frame); // starts at 1 !!!
  // keep only the begining of the prom
  int left_after = PROM_SEQ_L - nb_elem_in_first_frame;

  // rebuild the seg in a int32
  uint32_t seqAtPos;
  seqAtPos = (seq_.seq[frame] & ((1 << nb_elem_in_first_frame) -1)) << left_after; // mask and left shift
  for( ; left_after >=  FRAME_SIZE; left_after -= FRAME_SIZE)
    seqAtPos += seq_.seq[++frame] << (left_after - FRAME_SIZE); // just read and left shift
  seqAtPos += seq_.seq[++frame] >> (FRAME_SIZE - left_after); // right shift

  // hamming weight of the XOR
  uint32_t diff = PROM_SEQ ^ seqAtPos;
  int dist = 0;
  for( ; diff; ++dist) diff &= diff - 1;

  return dist;
}

int Dna::terminator_at(int pos) {
  // Only works for eSIZE <= FRAME_SIZE + 1
  // TODO modulo et sortie de range...
  int eSIZE = 4;
  int eOFFSET = 10;
  int frame;
  int nb_elem_in_first_frame;
  int left_over;

  uint32_t seqAtPos;
  frame = FRAME(pos);
  nb_elem_in_first_frame = FRAME_SIZE - INDEXFL(pos, frame); // starts at 1 !!!
  left_over = eSIZE - nb_elem_in_first_frame;
  if (left_over > 0){
    seqAtPos = (seq_.seq[frame] >> nb_elem_in_first_frame) & ((1 << nb_elem_in_first_frame) - 1)
      << left_over; // shift right and mask and left shift
    seqAtPos+= seq_.seq[++frame] >> (FRAME_SIZE - left_over);
  } else {
    seqAtPos = (seq_.seq[frame] >> (-left_over)) 
      & ((1 << eSIZE) - 1); // shift right and mask
  }

  uint32_t seqAtOffset;
  frame = FRAME(pos+eOFFSET);
  nb_elem_in_first_frame = FRAME_SIZE - INDEXFL(pos+eOFFSET, frame); // starts at 1 !!!
  left_over = eSIZE - nb_elem_in_first_frame;
  if (left_over > 0){
    seqAtOffset = (seq_.seq[frame] >> nb_elem_in_first_frame) & ((1 << nb_elem_in_first_frame) - 1)
      << left_over; // shift right and mask and left shift
    seqAtOffset+= seq_.seq[++frame] >> (FRAME_SIZE - left_over);
  } else {
    seqAtOffset = (seq_.seq[frame] >> (-left_over)) 
      & ((1 << eSIZE) - 1); // shift right and mask
  }

  uint32_t diff = seqAtPos ^ seqAtOffset;
  int dist = 0;
  for( ; diff; ++dist) diff &= diff - 1;

  return dist;
}

bool Dna::shine_dal_start(int pos) {
  // TODO modulo et sortie de range...
  bool start = false;

  int frame = FRAME(pos);
  int index = INDEXFL(pos, frame);
  int t_frame, t_index;

  for (int k = 0; k < 9; k++) {
    t_index = index + (k >= 6 ? k + 4 : k);
    for(t_frame = frame; t_index > FRAME_SIZE ; ++t_frame, t_index -= FRAME_SIZE);

    // XOR then mask on last bit
    bool diff = 0b1 & (
      (seq_.seq[t_frame] >> (FRAME_SIZE - t_index - 1)) 
      ^
      SHINE_DAL_SEQ >> (SHINE_DAL_SEQ_L - k - 1)
    );

    if (diff) {
      start = false;
      break;
    } else {
      start = true;
    }
  }

  return start;
}

bool Dna::protein_stop(int pos) {
  bool is_protein;

  int frame = FRAME(pos);
  int index = INDEXFL(pos, frame);
  int t_frame, t_index;

  for (int k = 0; k < 3; k++) {
    t_index = index + k;
    for(t_frame = frame; t_index > FRAME_SIZE ; ++t_frame, t_index -= FRAME_SIZE);
    
    // XOR then mask on last bit
    bool diff = 0b1 & (
      (seq_.seq[t_frame] >> (FRAME_SIZE - t_index - 1)) 
      ^
      PROTEIN_END >> (PROTEIN_END_L - k - 1)
    );

    if (diff) {
      is_protein = false;
      break;
    } else {
      is_protein = true;
    }
  }

  return is_protein;
}

int Dna::codon_at(int pos) {
  int value = 0;

  int frame = FRAME(pos);
  int index = INDEXFL(pos, frame);
  int t_frame, t_index;

  for (int i = 0; i < 3; i++) {
    t_index = index + i;
    for(t_frame = frame; t_index > FRAME_SIZE ; ++t_frame, t_index -= FRAME_SIZE);

    bool isOne = 0b1 & (seq_.seq[t_frame] >> (FRAME_SIZE - t_index - 1));
    if (isOne)
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}