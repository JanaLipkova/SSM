#ifndef ICOS_RNG_UNIFORM_MT_H
#define ICOS_RNG_UNIFORM_MT_H

typedef struct{
    unsigned int matrix_a;
    unsigned int mask_b;
    unsigned int mask_c;
    unsigned int seed;
} mt_struct_stripped;


#define N_GPU_RNG_BLOCKS				32
#define N_GPU_RNG_THREADS_PER_BLOCK		128

#define N_GPU_RNG_THREADS_TOTAL		N_GPU_RNG_BLOCKS * N_GPU_RNG_THREADS_PER_BLOCK
#define N_MT_STATES					MT_RNG_COUNT * MT_NN

#define GPU_RNG_INITIALSTATE		777

#define          MT_MM 9
#define          MT_NN 19
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0xFFFFFFFEU
#define       MT_LMASK 0x1U
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18


#define		MAX_N_RAND 4096

#endif