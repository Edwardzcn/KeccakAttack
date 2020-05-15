#pragma once
#pragma once
//#include <stdlib.h>
//#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
//#define reducedRound 3

typedef unsigned char UINT8;
typedef unsigned long long UINT64;
typedef UINT64 tKeccakLane;

#define DEBUG
#define KeccakReferences
#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define KeccakP1600_stateSizeInBytes    200
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)

static tKeccakLane KeccakRoundConstants[maxNrRounds];
static unsigned int KeccakRhoOffsets[nrLanes];

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
void score(UINT64* hash, int outputLen, int reducedRound);

static void theta(tKeccakLane* A);
static void rho(tKeccakLane* A);
static void pi(tKeccakLane* A);
static void chi(tKeccakLane* A);
static void iota(tKeccakLane* A, unsigned int indexRound);

UINT64 stringtoull(string str_64);

void test_preimage(string input_path);
