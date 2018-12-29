//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_.length){ 
  memcpy(seq_.seq, clone.seq_.seq, seq_.size);
}

Dna::Dna(int length, Threefry::Gen& rng) : seq_(length) {

  uint32_t max = 0 - 1;
  for (int32_t i = 0; i < seq_.nbElem - REPETITION_MARGIN; i++) {
    seq_.seq[i] = rng.random(max);
  }
  for (int32_t i = seq_.nbElem - REPETITION_MARGIN; i < seq_.nbElem; i++) {
    seq_.seq[i] = 0;
    for(int32_t j = 0; j < FRAME_SIZE; ++j){
      int pos = i * FRAME_SIZE + j;
      if(pos < length)
      {
        // random generate
        seq_.seq[i] += LEFT(rng.random(NB_BASE), (FRAME_SIZE - j - 1));
      }
      else
      {
        // repeat
        int frame = FRAME(pos - length);
        int p_in_frame_from_r = INDEXFR(pos - length, frame);
        seq_.seq[i] += 
          LEFT(
            AT(seq_.seq[frame], p_in_frame_from_r),
            (FRAME_SIZE - j - 1));
      }
    }
  }
  // Old seed generator (slower) :
  /*
  for (int32_t i = 0; i < seq_.nbElem; i++) {
    seq_.seq[i] = 0;
    for(int32_t j = 0; j < FRAME_SIZE; ++j){
      int pos = i * FRAME_SIZE + j;
      if(pos < length)
      {
        // random generate
        seq_.seq[i] += rng.random(NB_BASE) << (FRAME_SIZE - j - 1);
      }
      else
      {
        // repeat
        int frame = FRAME(pos - length);
        int p_in_frame_from_r = INDEXFR(pos - length, frame);
        seq_.seq[i] += 
          (((seq_.seq[frame] >> p_in_frame_from_r)
            & 0b1)
          << (FRAME_SIZE - j - 1));
      }
    }
  }
  */
}

Dna::Dna(uint32_t* genome, int length) : seq_(length) {
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
  uint32_t mask = LEFT(1, p_in_frame_from_r);

  if(c) seq_.seq[frame] |= mask; // set 1
  else seq_.seq[frame] &= ~mask; // set 0
}

void Dna::do_switch(int pos) {
  int frame = FRAME(pos);
  int p_in_frame_from_r = INDEXFR(pos, frame);
  uint32_t mask = LEFT(1, p_in_frame_from_r);

  seq_.seq[frame] ^= mask;
}

int Dna::promoter_at(int pos) {
  int frame = FRAME(pos);
  int nb_elem_in_first_frame = FRAME_SIZE - INDEXFL(pos, frame); // starts at 1 !!!
  // keep only the begining of the prom
  int left_after = PROM_SEQ_L - nb_elem_in_first_frame;

  // rebuild the seg in a int32
  uint32_t seqAtPos;
  if(left_after > 0){
    seqAtPos = LEFT(LAST(seq_.seq[frame], nb_elem_in_first_frame), left_after); // mask and left shift
    /*for( ; left_after >  FRAME_SIZE; left_after -= FRAME_SIZE)
      seqAtPos += LEFT(seq_.seq[++frame], left_after - FRAME_SIZE); // just read and left shift*/
    seqAtPos += FIRST(seq_.seq[++frame], left_after, FRAME_SIZE); // right shift
  } else {
    seqAtPos = LAST(RIGHT(seq_.seq[frame], -left_after), PROM_SEQ_L);
  }

  // hamming weight of the XOR
  uint32_t diff = PROM_SEQ ^ seqAtPos;
  int dist = 0;
  COUNT(diff, dist);

  return dist;
}

int Dna::terminator_at(int pos) {
  // Only works for eSIZE <= FRAME_SIZE + 1
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
    seqAtPos = LEFT(LAST(seq_.seq[frame], nb_elem_in_first_frame), left_over); // mask and left shift
    seqAtPos+= FIRST(seq_.seq[++frame], left_over, FRAME_SIZE);
  } else {
    seqAtPos = LAST(RIGHT(seq_.seq[frame], -left_over), eSIZE); // shift right and mask
  }

  int index2 = pos + eOFFSET - eSIZE + 1;
  frame = FRAME(index2);
  int p_in_frame_from_r = INDEXFR(index2, frame);
  uint32_t seqAtOffset = 0;
  for(int i = 0; i < eSIZE; ++i){
    seqAtOffset += LEFT(AT(seq_.seq[frame], p_in_frame_from_r), i);
    p_in_frame_from_r--;
    if(p_in_frame_from_r < 0)
      frame++;
  }

  uint32_t diff = seqAtPos ^ seqAtOffset;
  int dist = 0;
  COUNT(diff, dist);

  return dist;
}

bool Dna::shine_dal_start(int pos) {
  bool start = false;

  int frame = FRAME(pos);
  int index = INDEXFL(pos, frame);
  int t_frame, t_index;

  for (int k = 0; k < 9; k++) {
    t_index = index + (k >= 6 ? k + 4 : k);
    t_frame = frame;

    if(t_index >= FRAME_SIZE){
      t_frame++;
      t_index -= FRAME_SIZE;
    }

    // XOR
    bool diff = AT(seq_.seq[t_frame],FRAME_SIZE - t_index - 1)
      ^ AT(SHINE_DAL_SEQ, SHINE_DAL_SEQ_L - k - 1);

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
    t_frame = frame;

    if(t_index >= FRAME_SIZE){
      t_frame++;
      t_index -= FRAME_SIZE;
    }
    
    // XOR then mask on last bit
    bool diff = AT(seq_.seq[t_frame], FRAME_SIZE - t_index - 1)
      ^ AT(PROTEIN_END, PROTEIN_END_L - k - 1);

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
    t_frame = frame;

    if(t_index >= FRAME_SIZE){
      t_frame++;
      t_index -= FRAME_SIZE;
    }

    bool isOne = AT(seq_.seq[t_frame], FRAME_SIZE - t_index - 1);
    if (isOne)
      value += LEFT(1, CODON_SIZE - i - 1);
  }

  return value;
}