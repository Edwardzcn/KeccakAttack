#include "operation.h"

/******* 5 operation functions ************/

extern tKeccakLane KeccakRoundConstants[maxNrRounds];
extern unsigned int KeccakRhoOffsets[nrLanes];


#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)

void theta(tKeccakLane* A)
{
	unsigned int x, y;
	tKeccakLane C[5], D[5];

	for (x = 0; x < 5; x++) {
		C[x] = 0;
		for (y = 0; y < 5; y++)
			C[x] ^= A[index(x, y)];
	}
	for (x = 0; x < 5; x++)
		D[x] = ROL64(C[(x + 1) % 5], 1) ^ C[(x + 4) % 5];
	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			A[index(x, y)] ^= D[x];
}

void rho(tKeccakLane* A)
{
	unsigned int x, y;

	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		A[index(x, y)] = ROL64(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

void pi(tKeccakLane* A)
{
	unsigned int x, y;
	tKeccakLane tempA[25];

	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		tempA[index(x, y)] = A[index(x, y)];
	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		A[index(0 * x + 1 * y, 2 * x + 3 * y)] = tempA[index(x, y)];
}

void chi(tKeccakLane* A)
{
	unsigned int x, y;
	tKeccakLane C[5];

	for (y = 0; y < 5; y++) {
		for (x = 0; x < 5; x++)
			C[x] = A[index(x, y)] ^ ((~A[index(x + 1, y)]) & A[index(x + 2, y)]);
		for (x = 0; x < 5; x++)
			A[index(x, y)] = C[x];
	}
}

void chi_inverse(tKeccakLane* A)
{
	unsigned int x, y;
	tKeccakLane C[5];

	for (y = 0; y < 5; ++y) {
		for (x = 0; x < 5; ++x) {
			C[x] = A[index(x, y)] ^ ((~A[index(x + 1, y)]) & (A[index(x + 2, y)] ^ ((~A[index(x + 3, y)]) & A[index(x + 4, y)])));
		}
		for (x = 0; x < 5; ++x)
			A[index(x, y)] = C[x];
	}
}

void iota(tKeccakLane* A, unsigned int indexRound)
{
	A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}
void iota_inverse(tKeccakLane* A, unsigned int indexRound) {
	A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}


#include <iostream>

void iota_inverse(vector<vec_GF2>& A, unsigned int indexRound, tKeccakLane* rc) {
	tKeccakLane value = rc[indexRound];
	for (int i = 0; i < 64; ++i) {
		if (value & 1LL << i) {
			cout << "# i=" << i << endl;
			A[0][i] = A[0][i] + 1;
		}
	}

}
/******* 5 operation functions ************/
