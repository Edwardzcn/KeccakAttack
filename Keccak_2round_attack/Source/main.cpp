#include "attack.h"
#include "verify.h"
#include <iostream>
#include <sstream>


using namespace std;

int main() {

	string 	input_path = string(".\\Data\\input_hash_value.txt");
	string  output_path = string(".\\Data\\output_preimage.txt");
	get_preimage(input_path, output_path);
	test_preimage(output_path);
}