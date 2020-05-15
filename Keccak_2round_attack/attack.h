#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include <string>
#include <vector>

void output_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length, ofstream& of);
void output_print_lane(vector<vec_GF2> &g, char vec_name, int bit_length, int lane_length, ofstream& of);
void debug_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length);
void debug_print_lane(vector<vec_GF2> &g, char vec_name, int bit_length, int lane_length);

string vecGF2tostring(vec_GF2 g, int bit_length);
//UINT64 stringtoull(string str_64);

void get_preimage(string input_path, string output_path);


vector<vec_GF2> get192zero(int bitsize, vec_GF2* hh);
vector<vec_GF2> get_test_input(int bitsize, vec_GF2* hh);
vector<vec_GF2> get_test_input_2(int bitsize, vec_GF2* hh);
