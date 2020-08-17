#include "attack.h"
#include "verify.h"
#include <iostream>
#include <sstream>


using namespace std;

int main() {

	string 	input_path = string(".\\Data\\input_hash_value.txt");
	string  output_path = string(".\\Data\\output_preimage.txt");
	get_preimage(input_path, output_path);
	//test_preimage_ull(output_path);
	independent_bits_engine<default_random_engine, 64, unsigned long long int> engine;
	engine.discard(30000000);
	rep(i, 0, 10000000) {
		test_preimage_random(0,engine);
	}
}