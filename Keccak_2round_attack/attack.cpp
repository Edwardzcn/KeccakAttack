#include "attack.h"

NTL_CLIENT

using namespace std;

void debug_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length) {
	cout << "***************** MESSAGE "<<vec_name<<" ****************\n";
	rep(i, 0, lane_length) {
		cout << vec_name << "[" << i << "]: ";
		rep(j, 0, bit_length) {
			cout << g[i][j];
		}
		cout << endl;
	}
	cout << "\n****************************************\n\n";
}

void debug_print_lane(vector<vec_GF2> &g, char vec_name, int bit_length, int lane_length) {
	cout << "***************** MESSAGE " << vec_name << " ****************\n";
	rep(i, 0, lane_length) {
		cout << vec_name << "[" << i << "]: ";
		rep(j, 0, bit_length) {
			cout << g[i][j];
		}
		cout << endl;
	}
	cout << "\n****************************************\n\n";
}

void output_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length, ofstream& of) {
	rep(i, 0, lane_length) {
		//of << vec_name << "[" << i << "]: ";
		rep(j, 0, bit_length) {
			of << g[i][j];
		}
		of << endl;
	}
}

void output_print_lane(vector<vec_GF2>& g, char vec_name, int bit_length, int lane_length, ofstream& of) {
	rep(i, 0, lane_length) {
		//of << vec_name << "[" << i << "]: ";
		rep(j, 0, bit_length) {
			of << g[i][j];
		}
		of << endl;
	}
}

string vecGF2tostring(vec_GF2 g, int bit_length) {
	stringstream pipe;
	string transtr;
	if (bit_length == 64) {
		pipe.clear();
		rep(i, 0, bit_length) {
			pipe << g[i];
		}
		pipe >> transtr;
		return transtr;
	}
	else {
		return "";
	}
}

//UINT64 stringtoull(string str_64) {
//	// 2进制字符串 转UINT64
//	return std::strtoull(str_64.c_str(), NULL, 2);
//}

void get_preimage(string input_path, string output_path)
{

	/*****************************************************/
	//Opening the file which contains the hash value.
	/*****************************************************/

	char ch;

	int i, j;

	vec_GF2 h[8];
	for (i = 0; i < 8; ++i) {
		h[i].SetLength(64);
		//		random(h[i], 64);	//change this
	}

	cout << "The filename to be opened is : " << input_path;


	ifstream infile;
	infile.open(input_path);
	if (!infile) {
		printf("Cannot open file \n");
		exit(0);
	}

	i = 0;
	infile.get(ch);
	while (i < 512) {
		h[i / 64][i % 64] = ch - '0';
		++i;
		infile.get(ch);
	}

	infile.close();

	/*****************************************************/
	//End of reading file
	/*****************************************************/

	h[0][0] = h[0][0] + 1; // iota inverse

	vec_GF2 hh[9];
	for (i = 0; i < 9; ++i)
		hh[i].SetLength(64);

	for (i = 0; i < 64; ++i) // chi inverse
	{
		hh[8][i] = 1;
		hh[7][i] = h[7][i];
		hh[6][i] = (h[6][i] + h[7][i] + 1);
		hh[5][i] = (h[5][i] + (h[6][i] + h[7][i]) * h[7][i]);

		hh[0][i] = (h[0][i] + (h[1][i] + 1) * (h[2][i] + (h[3][i] + 1) * h[4][i]));
		hh[1][i] = (h[1][i] + (h[2][i] + 1) * (h[3][i] + (h[4][i] + 1) * h[0][i]));
		hh[2][i] = (h[2][i] + (h[3][i] + 1) * (h[4][i] + (h[0][i] + 1) * h[1][i]));
		hh[3][i] = (h[3][i] + (h[4][i] + 1) * (h[0][i] + (h[1][i] + 1) * h[2][i]));
		hh[4][i] = (h[4][i] + (h[0][i] + 1) * (h[1][i] + (h[2][i] + 1) * h[3][i]));
	}
	// 作者为了省事  没有再另创建数组保存/iota^-1的运算结果
	h[0][0] = h[0][0] + 1; //to get back the hash value

	vector<vec_GF2> m_d, m_e;
	//int m_d_bitsize = 1088;
	//m_d = get192zero(1088,hh);
	m_d = get_test_input(1088, hh);
	m_e = get_test_input_2(1088, hh);


	cout << "****************HASH VALUE**************\n";

	// change hash size
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 64; ++j) {
			cout << h[i][j];
		}
		cout << endl;
	}
	cout << "****************************************\n\n";

	cout << "The preimage is stored in the file named \"preimage\".\n";
	ofstream outfile;
	outfile.open(output_path, ios::out);
	if (!outfile) {
		cout << "Cannot open file \n";
		exit(1);
	}


	debug_print_lane(m_d, 'd', 64, 17);
	output_print_lane(m_d, 'd', 64, 17, outfile);
	debug_print_lane(m_e, 'e', 64, 17);
	output_print_lane(m_e, 'e', 64,17, outfile);


	


	outfile.close();
}

vector<vec_GF2> get192zero(int bitsize,vec_GF2* hh)
{
	int lanesize = bitsize / 64;
	vector<vec_GF2> e;
	for (int i = 0; i < 17; i++)
	{
		e.push_back(vec_GF2());
		e[i].SetLength(64);
	}


	for (int i = 0; i < 64; ++i) {
		e[0][i] = hh[0][i];  // e_0 e_5  e_10  e_15 
		e[5][i] = hh[0][i];
		e[10][i] = hh[0][i];
		e[15][i] = hh[0][i];
		//e[1][i] = e[6][i] = e[11][i] = e[16][i] = 0;  // e_1  e_6 e_11 e_16   
		//e[2][i] = e[7][i] = e[12][i] = 0;
		//e[3][i] = e[8][i] = e[13][i] = 0;
	}

	//for padding

	e[11][0] = e[11][63] = e[16][0] = e[16][63] = 1;

	/*****************************************************/
	//End
	/*****************************************************/

	return e;
}



vector<vec_GF2> get_test_input(int bitsize, vec_GF2* hh)
{
	int lanesize = bitsize / 64;
	vector<vec_GF2> d;
	for (int i = 0; i < 17; i++)
	{
		d.push_back(vec_GF2());
		d[i].SetLength(64);
	}


	//for (int i = 0; i < 64; ++i) {

	//}

	//for padding
	// 构造思路  
	// 一次错误构造：未考虑 h' 偏移
	//d[4][64-27] = d[9][64-27] = d[1][64-45] = d[16][64-45] = 1;

	d[4][(64 - 27 + 43 + 64) % 64] = d[9][(64 - 27 +43 + 64) % 64] = d[1][64 - 45] = d[16][64 - 45] = 1;
		
	/*****************************************************/
	//End
	/*****************************************************/

	return d;
}

vector<vec_GF2> get_test_input_2(int bitsize, vec_GF2* hh)
{
	int lanesize = bitsize / 64;
	vector<vec_GF2> e;
	for (int i = 0; i < 17; i++)
	{
		e.push_back(vec_GF2());
		e[i].SetLength(64);
	}


	//for (int i = 0; i < 64; ++i) {

	//}

	//for padding

	//e[1][0] = e[1][63] = e[11][0] = e[9][7] = 1;
	//e[3][1] = e[10][44] = e[13][44] = e[7][7] = e[8][0] = 1;


	//e[1][0] = e[1][64-63] = e[11][0] = e[9][64-7] = 1;
	//e[3][64-1] = e[10][64-44] = e[13][64-44] = e[7][64-7] = e[8][0] = 1;

	//e[1][0] = e[1][63] = 1;
	//e[10][20] = e[6][57] = e[9][57] = 1;
	//e[11][0] = e[13][20] = 1;


	//// x=0
	//e[5][0] = e[10][20] = e[15][21] = 1;
	//// x=1
	//e[1][44] = e[6][0] = e[6][14] = e[6][44] = 1;
	//e[11][0] = e[11][63] = 1;
	//// x=2
	//// pass
	//// x=3
	//e[8][0] = e[13][20] = e[13][21] = 1;
	//// x=4
	//e[9][14] = 1;


	// x=0
	e[5][0] = e[10][20] = e[15][43] = 1;
	// x=1
	e[6][36] = e[6][0] = e[6][20] = 1;
	e[1][20] = 1;
	e[11][0] = e[11][63] = 1;
	// x=2
	// pass
	// x=3
	e[8][0] = e[13][20] = e[13][43] = 1;
	// x=4
	e[9][36] = 1;


	/*****************************************************/
	//End
	/*****************************************************/

	return e;
}