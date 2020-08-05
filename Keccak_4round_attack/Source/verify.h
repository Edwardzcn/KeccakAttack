#pragma once
#include <memory.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <cstdlib>
#include "common.h"
#include "operation.h"

using namespace std;

typedef unsigned char UINT8;
typedef unsigned long long UINT64;
typedef UINT64 tKeccakLane;





void KeccakP1600_InitializeRoundConstants(void);
void KeccakP1600_InitializeRhoOffsets(void);
static int LFSR86540(UINT8* LFSR);
void KeccakP1600_InitializeRoundConstants(void);
void KeccakP1600_InitializeRhoOffsets(void);

void KeccakP1600Round(tKeccakLane* state, unsigned int indexRound);
void KeccakF1600Permute(tKeccakLane* state, unsigned int nRounds);

void Keccak(int reducedRound, unsigned int rate, unsigned int capacity, const UINT64* input,
	unsigned int inputBitLen, UINT64* output,
	unsigned int outputBitLen, int level);

void print1600state(tKeccakLane* state);
void printHash(UINT64 message[], int inputLen, UINT64* hash, int outputLen);

extern UINT64 stringtoull(string str_64);
void score(UINT64* hash, int outputLen, int reducedRound);
int score_return_len(UINT64* hash, int outputLen, int reducedRound);
void test_preimage(string input_path);
void test_preimage_ull(string input_path);
void test_preimage_random(int pid, independent_bits_engine<default_random_engine, 64, unsigned long long int>& input_engine);



