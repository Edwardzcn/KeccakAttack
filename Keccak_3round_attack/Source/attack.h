#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include <string>
#include <vector>
#include "operation.h"





#define INF 0x3f3f3f3f
#define HASH_BIT_SIZE 256


extern void GetRoundConstants(tKeccakLane* rc);
extern void print1600state(tKeccakLane* state);

void output_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length, ofstream& of);
void output_print_lane(vector<vec_GF2> &g, char vec_name, int bit_length, int lane_length, ofstream& of);
void debug_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length);
void debug_print_lane(vector<vec_GF2> &g, char vec_name, int bit_length, int lane_length);

string vecGF2tostring(vec_GF2 g, int bit_length);
UINT64 stringtoull(string str_64);

void get_preimage(string input_path, string output_path);


vector<vec_GF2> get192zero(int bitsize, vector<vec_GF2>& hh);
vector<vec_GF2> get_test_input(int bitsize, vector<vec_GF2>& hh);
vector<vec_GF2> get_test_input_2(int bitsize, vector<vec_GF2>& hh);

const int maxn = 400;


void gaussian_init(int a[maxn][maxn],int x[maxn],int &equ,int &var);
int gaussian_solve(int a[maxn][maxn], int x[maxn],int RC[5][maxn],int free_x[maxn],int &equ,int &var,int &free_num);


