#include "Algorithms.h"
#include "Algorithms.cuh"

#include "ExpManager.h"
#include "ThreefryGPU.h"
#include "GPUDna.cuh"

#include <cstdint>
#include <stdio.h>
#include <unistd.h>

#include <iostream>

#include<cuda.h>
#include<cuda_profiler_api.h>

using namespace std;

#define DEBUG 1
// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n",
                cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

#define SELECTION_SCOPE_X 3
#define SELECTION_SCOPE_Y 3
#define HALF_SCOPE_X 1
#define HALF_SCOPE_Y 1
#define NEIGHBORHOOD_SIZE 9


constexpr int32_t PROMOTER_ARRAY_SIZE = 10000;
uint32_t* gpu_sequences;
uint32_t* gpu_prev_sequences;
double* gpu_fitness;

/**
 * Function to transfer data from CPU to GPU
 *
 * @param exp_m
 * @param first_gen
 */
 void transfer_in(ExpManager* exp_m, bool first_gen) {
     exp_m->rng_->initDevice();

     // Malloc
     checkCuda(cudaMalloc((void**) &gpu_counters,
                          exp_m->rng_->counters().size() *
                          sizeof(unsigned long long)));
    checkCuda(cudaMalloc((void**) &gpu_prev_sequences,
                        exp_m->nb_indivs_ * 
                        exp_m->internal_organisms_[0]->dna_->seq_.size));
    
    checkCuda(cudaMalloc((void**) &gpu_sequences, sizeof(gpu_prev_sequences)));
    
    checkCuda(cudaMalloc((void**) &gpu_fitness, exp_m->nb_indivs_ * sizeof(double)));
    
    // Mem cpy
    checkCuda(cudaMemcpy(gpu_counters, exp_m->rng_->counters().data(),
                        exp_m->rng_->counters().size() *
                        sizeof(unsigned long long), cudaMemcpyHostToDevice));                    
    for(int i = 0; i < exp_m->nb_indivs_; ++i) {
        checkCuda(cudaMemcpy(&gpu_prev_sequences[i], exp_m->internal_organisms_[i]->dna_->seq_.seq,
                            exp_m->internal_organisms_[0]->dna_->seq_.size, 
                            cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(&gpu_fitness[i], &(exp_m->internal_organisms_[i]->fitness),
                            sizeof(double), 
                            cudaMemcpyHostToDevice));
        //printf("%f - ", exp_m->internal_organisms_[i]->fitness);
    }

     // TO COMPLETE
 }

 void transfer_out(ExpManager* exp_m) {
     // TODO
 }

 void clean(ExpManager* exp_m) {
    checkCuda(cudaFree(gpu_counters));
    checkCuda(cudaFree(gpu_sequences));
    checkCuda(cudaFree(gpu_prev_sequences));
    // TO COMPLETE
}

__device__ int32_t Threefry::Device::roulette_random(double* probs, int32_t nb_elts)
{
    double pick_one = 0.0;

    while (pick_one == 0.0)
    {
        pick_one = randomDouble();
    }

    int32_t found_org = 0;

    pick_one -= probs[0];
    while (pick_one > 0)
    {
        assert(found_org<nb_elts-1);

        pick_one -= probs[++found_org];
    }
    return found_org;
}


__constant__ double cof[6] = {  76.18009172947146,
                                -86.50532032941677,
                                24.01409824083091,
                                -1.231739572450155,
                                0.1208650973866179e-2,
                                -0.5395239384953e-5 };



// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt).
// When the argument z is an integer, the gamma function is just the familiar factorial
// function, but offset by one, n! = gamma(n + 1).
__device__ static double gammln(double X)
{
    double x, y, tmp, ser;

    y = x = X;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;

    for (int8_t j = 0 ; j <= 5 ; j++)
    {
        ser += cof[j] / ++y;
    }

    return -tmp + log(2.5066282746310005 * ser / x);
}


__device__ 
int32_t Threefry::Device::binomial_random(int32_t nb_drawings, double prob)
{
    int32_t nb_success;

    // The binomial distribution is invariant under changing
    // ProbSuccess to 1-ProbSuccess, if we also change the answer to
    // NbTrials minus itself; we ll remember to do this below.
    double p;
    if (prob <= 0.5) p = prob;
    else p = 1.0 - prob;

    // mean of the deviate to be produced
    double mean = nb_drawings * p;


    if (nb_drawings < 25)
        // Use the direct method while NbTrials is not too large.
        // This can require up to 25 calls to the uniform random.
    {
        nb_success = 0;
        for (int32_t j = 1 ; j <= nb_drawings ; j++)
        {
            if (randomDouble() < p) nb_success++;
        }
    }
    else if (mean < 1.0)
        // If fewer than one event is expected out of 25 or more trials,
        // then the distribution is quite accurately Poisson. Use direct Poisson method.
    {
        double g = exp(-mean);
        double t = 1.0;
        int32_t j;
        for (j = 0; j <= nb_drawings ; j++)
        {
            t = t * randomDouble();
            if (t < g) break;
        }

        if (j <= nb_drawings) nb_success = j;
        else nb_success = nb_drawings;
    }

    else
        // Use the rejection method.
    {
        double en     = nb_drawings;
        double oldg   = gammln(en + 1.0);
        double pc     = 1.0 - p;
        double plog   = log(p);
        double pclog  = log(pc);

        // rejection method with a Lorentzian comparison function.
        double sq = sqrt(2.0 * mean * pc);
        double angle, y, em, t;
        do
        {
            do
            {
                angle = M_PI * randomDouble();
                y = tan(angle);
                em = sq*y + mean;
            } while (em < 0.0 || em >= (en + 1.0)); // Reject.

            em = floor(em); // Trick for integer-valued distribution.
            t = 1.2 * sq * (1.0 + y*y)
                * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);

        } while (randomDouble() > t); // Reject. This happens about 1.5 times per deviate, on average.

        nb_success = (int32_t) rint(em);
    }


    // Undo the symmetry transformation.
    if (p != prob) nb_success = nb_drawings - nb_success;

    return nb_success;
}

__device__ static int mod(int a, int b)
{

    assert(b > 0);

    while (a < 0)  a += b;
    while (a >= b) a -= b;

    return a;
}

__global__ void selection(unsigned long long* gpu_counters, double* fitness, uint32_t* seqs, uint32_t* next,
    int nb_indiv, int grid_w, int grid_h, int size_seq){
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int indiv_id = x * grid_h + y;

    if(x < grid_w && y < grid_h){
        double probs[NEIGHBORHOOD_SIZE];
        double tot_local_probs = 0;
        int count = 0;
        int cur_x;
        int cur_y;
        for(int i = - HALF_SCOPE_X; i  < SELECTION_SCOPE_X - HALF_SCOPE_X; ++i){
            for(int j = - HALF_SCOPE_Y; j < SELECTION_SCOPE_Y - HALF_SCOPE_Y; ++j){
                cur_x = mod(x+i, grid_w);
                cur_y = mod(y+j, grid_h);
    
                probs[count] = fitness[cur_x * grid_h + cur_y];
                tot_local_probs +=  probs[count];
    
                ++count;
            }
        }
        for(int i = 0 ; i < NEIGHBORHOOD_SIZE ; ++i) {
            probs[i] = probs[i]/tot_local_probs;
        }
    
        Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::REPROD,nb_indiv);
        int found_org = rng.roulette_random(probs, NEIGHBORHOOD_SIZE);
    
        int x_offset = (found_org / SELECTION_SCOPE_X) - HALF_SCOPE_X;
        int y_offset = (found_org % SELECTION_SCOPE_Y) - HALF_SCOPE_Y;
        cur_x = mod(x+x_offset, grid_w);
        cur_y = mod(y+y_offset, grid_h);
        int next_indiv_id = cur_x * grid_h + cur_y;
    
        for(int i = 0 ; i < size_seq ; ++i) {
            next[indiv_id * size_seq + i] = seqs[next_indiv_id * size_seq + i];
        }
    }
}
// seq_length is its true length in number of bool
__global__ void do_mutation(unsigned long long* gpu_counters, uint32_t* seqs,
    int nb_indiv, int grid_w, int grid_h, int seq_length, int size_seq, 
    double mutation_rate) {

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int indiv_id = x * grid_h + y;

    if(x < grid_w && y < grid_h){
        Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::MUTATION,nb_indiv);
        int nb_mut = rng.binomial_random(seq_length, mutation_rate);
        if(nb_mut > 0) {
            int nb_swi = nb_mut;
    
            while(nb_mut > 0) {
                int rand_val = rng.random(nb_mut);
                if(rand_val < nb_swi){
                    --nb_swi;
                    int pos = rng.random(seq_length);
                    dna_gpu_do_switch(seqs, indiv_id, size_seq, pos);
                }
                --nb_mut;
            }
        }
    }
}

/**
 * Run a step on the GPU
 * @param nb_indiv
 * @param seq_length number of bool in the seq of one indiv
 * @param w_max
 * @param selection_pressure
 * @param grid_width
 * @param grid_height
 * @param mutation_rate
 */
void run_a_step_on_GPU(int nb_indiv, int size_seq, int seq_length, double w_max, double selection_pressure, int grid_width, int grid_height, double mutation_rate) {
    dim3 DimGridOrganism(ceil(grid_width/16.),ceil(grid_height/16.),1);
    dim3 DimBlockOrganism(16,16,1);

    selection<<<DimGridOrganism, DimBlockOrganism>>>(gpu_counters, gpu_fitness, gpu_prev_sequences, gpu_sequences, 
        nb_indiv, grid_width, grid_height, size_seq);
    checkCuda(cudaGetLastError());

    do_mutation<<<DimGridOrganism, DimBlockOrganism>>>(gpu_counters, gpu_sequences, 
        nb_indiv, grid_width, grid_height, seq_length, size_seq, mutation_rate);
    checkCuda(cudaGetLastError());

    apply_next_gen();
}

/**
 * Reallocate some data structures if needed
 * @param nb_indiv
 */
void apply_next_gen() {
    uint32_t* temp = gpu_prev_sequences;
    gpu_prev_sequences = gpu_sequences;
    gpu_sequences = temp;
}

/**
PRNG usage:
 * For selection
        Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::REPROD,nb_indiv);
        int found_org = rng.roulette_random(probs, NEIGHBORHOOD_SIZE);
 * For mutation:
      Threefry::Device rng(gpu_counters,indiv_id,Threefry::Phase::MUTATION,nb_indiv);
      rng.binomial_random(prev_gen_size, mutation_r);
      rng.random( number );
 **/