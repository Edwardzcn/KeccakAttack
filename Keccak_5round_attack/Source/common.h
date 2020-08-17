#pragma once

#include <NTL/GF2.h>
#include <NTL/ZZ.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/vector.h>


#define rep(i,a,n) for(int i=a;i<n;++i)
#define per(i,a,n) for(int i=n-1;i>=a;--i)
#define fi first
#define se second

// #define DEBUG
#define ReducedRound 5
#define KeccakReferences
#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define KeccakP1600_stateSizeInBytes    200
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)
#define TARGET_LEN 22

typedef unsigned char UINT8;
typedef unsigned long long UINT64;
typedef UINT64 tKeccakLane;

NTL_CLIENT




//static tKeccakLane KeccakRoundConstants[maxNrRounds];
//static unsigned int KeccakRhoOffsets[nrLanes];

