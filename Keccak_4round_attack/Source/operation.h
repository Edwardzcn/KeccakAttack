#pragma once
#include "common.h"
#include <vector>


void theta(tKeccakLane* A);
void rho(tKeccakLane* A);
void pi(tKeccakLane* A);
void chi(tKeccakLane* A);
void chi_inverse(tKeccakLane* A);
void iota(tKeccakLane* A, unsigned int indexRound);
void iota_inverse(tKeccakLane* A, unsigned int indexRound);
void iota_inverse(vector<vec_GF2>& A, unsigned int indexRound, tKeccakLane* rc);


