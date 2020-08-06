#include "attack.h"
#include <windows.h>


NTL_CLIENT

using namespace std;

void debug_print_lane(vec_GF2* g, char vec_name, int bit_length, int lane_length) {
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

void debug_print_lane(vector<vec_GF2>& g, char vec_name, int bit_length, int lane_length) {
	cout << "***************** MESSAGE " << vec_name << " ****************\n";
	rep(i, 0, lane_length) {
		cout << vec_name << "[" << i << "]: ";
		rep(j, 0, bit_length) {
			cout << g[i][j];
		}
		cout << endl;
	}
	cout << "****************************************\n\n";
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

UINT64 stringtoull(string str_64) {
	// 2进制字符串 转UINT64
	return std::strtoull(str_64.c_str(), NULL, 2);
}

void get_preimage(string input_path, string output_path)
{

	/*****************************************************/
	//Opening the file which contains the hash value.
	/*****************************************************/

	char ch;
	int i, j;

	vector<vec_GF2> h;
	const int hash_lane_num = HASH_BIT_SIZE / 64 + 1;
	for (i = 0; i < hash_lane_num; ++i) {
		h.push_back(vec_GF2());
		h[i].SetLength(64);
	}


//#ifdef DEBUG
//	cout << "The filename to be opened is : " << input_path << endl;
//#endif // DEBUG

	//ifstream infile;
	//infile.open(input_path);
	//if (!infile) {
	//	printf("Cannot open file \n");
	//	exit(0);
	//}

	//i = 0;
	//infile.get(ch);
	//while (i < HASH_BIT_SIZE) {
	//	h[i / 64][i % 64] = ch - '0';
	//	++i;
	//	infile.get(ch);
	//}

	//infile.close();

	/*****************************************************/
	//End of reading file
	/*****************************************************/










	// 做RC的预处理
	tKeccakLane rc[25];
	GetRoundConstants(rc);





	vector<vec_GF2> hh(h);
//	vector<vec_GF2> hh;
//	for (i = 0; i < hash_lane_num; ++i) {
//		hh.push_back(vec_GF2());
//		hh[i].SetLength(64);
//		//		random(h[i], 64);	//change this  随机数赋值
//	}
//	for (i = 0; i < 64; ++i) // chi inverse
//	{
//		// get hh'
//		hh[0][i] = (hiota[0][i] + (hiota[1][i] + 1) * (hiota[2][i] + (hiota[3][i] + 1) * hiota[4][i]));
//		hh[1][i] = (hiota[1][i] + (hiota[2][i] + 1) * (hiota[3][i] + (hiota[4][i] + 1) * hiota[0][i]));
//		hh[2][i] = (hiota[2][i] + (hiota[3][i] + 1) * (hiota[4][i] + (hiota[0][i] + 1) * hiota[1][i]));
//		hh[3][i] = (hiota[3][i] + (hiota[4][i] + 1) * (hiota[0][i] + (hiota[1][i] + 1) * hiota[2][i]));
//		hh[4][i] = (hiota[4][i] + (hiota[0][i] + 1) * (hiota[1][i] + (hiota[2][i] + 1) * hiota[3][i]));
//	}
//
//#ifdef DEBUG
//	debug_print_lane(hh, 'c', 64, 5);
//#endif // DEBUG

	iota_inverse(hh, 2, rc);
#ifdef DEBUG
	cout << "After Create RC and final hiota (32bits) ************\n";
	debug_print_lane(hh, 'i', 64, 5);
#endif // DEBUG



	/*****************************************************/
	//Building a system of linear equation. Ax = b where x is a vector containing lane variables.
	/*****************************************************/

	vec_GF2 b, x;
	vec_GF2 RC0, RC1, RC2;
#define RANK_SIZE 384
	b.SetLength(RANK_SIZE);
	x.SetLength(RANK_SIZE);
	RC0.SetLength(64);
	RC1.SetLength(64);
	RC2.SetLength(64);
	for (int i = 0; i < 64; ++i) {
		if (rc[0] & 1LL << i) {
			RC0[i] = RC0[i] + 1;
		}
		if (rc[1] & 1LL << i) {
			RC1[i] = RC1[i] + 1;
		}
		if (rc[2] & 1LL << i) {
			RC2[i] = RC2[i] + 1;
		}
	}

	for (i = 0; i < 64; ++i) {

		// 注意求反当作常数
		//b[(0 * 64) + i] = hh[0][i] + RC[i] + 1;
		//b[(1 * 64) + i] = hh[1][(64 + i - 20) % 64] + RC[i];
		//b[(2 * 64) + i] = hh[2][(64 + i - 21) % 64] + 1;
		//b[(3 * 64) + i] = hh[3][(64 + i - 43) % 64];
		//b[(4 * 64) + i] = hh[4][(64 + i - 50) % 64] + RC[(64 + i - 1) % 64] + 1;
		//b[(5 * 64) + i] = 0;
		//b[(6 * 64) + i] = 0;
		b[0 * 64 + i] = 0;
		b[1 * 64 + i] = 0;
		b[2 * 64 + i] = RC0[i];
		b[3 * 64 + i] = 1;
		b[4 * 64 + i] = 0;
	}
	for (i = 0; i < 32; ++i) {
		b[320 + i] = RC1[i] + RC2[i] + 1;
		//b[352 + i] = RC0[(64 + i + 20) % 64] + RC1[(64 + i + 20) % 64] + 1;
		b[352 + i] = RC0[(64 + i + 63) % 64];
	}


	mat_GF2 A;
	// 构建 x[0] - x[6] 共7lanes矩阵
	//A.SetDims(448, 448); //猜测无错   设置A为384*384的矩阵
	A.SetDims(384, 384);
	clear(A);

	for (i = 0; i < 64; ++i) {
		//First equation
		//A[0 + i][(0 * 64) + (64 + i - 0) % 64] = 1; //x_0(0)
		//A[0 + i][(6 * 64) + ((64 + i - 43) % 64)] = 1; //x_6(43)
		//A[0 + i][(1 * 64) + (64 + i - 36) % 64] = 1; //x_1(36)
		//A[0 + i][(2 * 64) + ((64 + i - 4) % 64)] = 1; //x_2(4)
		//A[0 + i][(5 * 64) + ((64 + i - 7) % 64)] = 1; //x_5(7)
		//A[0 + i][(1 * 64) + ((64 + i - 37) % 64)] = 1; //x_1(37)
		//A[0 + i][(3 * 64) + ((64 + i - 42) % 64)] = 1; //x_3(42)
		A[0 + i][(0 * 64) + ((64 + i - 0) % 64)] = 1; //x_0(0)
		A[0 + i][(1 * 64) + ((64 + i - 0) % 64)] = 1; //x_1(0)
		A[0 + i][(2 * 64) + ((64 + i - 0) % 64)] = 1; //x_2(0)

		//Second equation
		//A[64 + i][(0 * 64) + ((64 + i - 0) % 64)] = 1; //x_0(0)
		//A[64 + i][(6 * 64) + ((64 + i - 43) % 64)] = 1; //x_6(43)
		//A[64 + i][(4 * 64) + ((64 + i - 62) % 64)] = 1; //x_4(62)
		//A[64 + i][(6 * 64) + ((64 + i - 44) % 64)] = 1; //x_6(44)
		//A[64 + i][(2 * 64) + ((64 + i - 4) % 64)] = 1; //x_2(4)
		A[64 + i][(3 * 64) + ((64 + i - 0) % 64)] = 1; //x_3(0)
		A[64 + i][(4 * 64) + ((64 + i - 0) % 64)] = 1; //x_4(0)
		A[64 + i][(5 * 64) + ((64 + i - 0) % 64)] = 1; //x_5(0)

		//Third equation
		//A[128 + i][(2 * 64) + ((64 + i - 3) % 64)] = 1; //x_2(3)
		//A[128 + i][(5 * 64) + ((64 + i - 6) % 64)] = 1; //x_5(6)
		//A[128 + i][(1 * 64) + ((64 + i - 36) % 64)] = 1; //x_1(36)
		//A[128 + i][(3 * 64) + ((64 + i - 41) % 64)] = 1; //x_3(41)
		//A[128 + i][(0 * 64) + ((64 + i - 1) % 64)] = 1; //x_0(1)
		//A[128 + i][(3 * 64) + ((64 + i - 42) % 64)] = 1; //x_3(42)
		//A[128 + i][(4 * 64) + ((64 + i - 63) % 64)] = 1; //x_4(63)
		A[128 + i][(0 * 64) + ((64 + i + 0) % 64)] = 1; //x_0(0)
		A[128 + i][(3 * 64) + ((64 + i + 2) % 64)] = 1; //x_3(2)
		A[128 + i][(5 * 64) + ((64 + i + 21) % 64)] = 1; //x_5(21)


		//Fourth equation
		//A[192 + i][(2 * 64) + ((64 + i - 3) % 64)] = 1; //x_2(3)
		//A[192 + i][(6 * 64) + ((64 + i - 43) % 64)] = 1; //x_6(43)
		//A[192 + i][(1 * 64) + ((64 + i - 37) % 64)] = 1; //x_1(37)
		A[192 + i][(1 * 64) + ((64 + i + 28) % 64)] = 1; //x_1(28)
		A[192 + i][(4 * 64) + ((64 + i + 58) % 64)] = 1; //x_4(58)


		//Fifth equation
		//A[256 + i][(0 * 64) + ((64 + i - 0) % 64)] = 1; //x_0(0)
		//A[256 + i][(3 * 64) + ((64 + i - 41) % 64)] = 1; //x_3(41)
		//A[256 + i][(4 * 64) + ((64 + i - 62) % 64)] = 1; //x_4(62)
		//A[256 + i][(0 * 64) + ((64 + i - 1) % 64)] = 1; //x_0(1)
		//A[256 + i][(6 * 64) + ((64 + i - 44) % 64)] = 1; //x_6(44)
		//A[256 + i][(2 * 64) + ((64 + i - 4) % 64)] = 1; //x_2(4)
		//A[256 + i][(4 * 64) + ((64 + i - 63) % 64)] = 1; //x_4(63)
		A[256 + i][(2 * 64) + ((64 + i + 61) % 64)] = 1; //x_2(61)
		A[256 + i][(5 * 64) + ((64 + i + 21) % 64)] = 1; //x_5(21)


		//Sixth equation
		//A[320 + i][(0 * 64) + ((64 + i - 0) % 64)] = 1; //x_0(0)
		//A[320 + i][(1 * 64) + ((64 + i - 0) % 64)] = 1; //x_1(0)
		//A[320 + i][(2 * 64) + ((64 + i - 0) % 64)] = 1; //x_2(0)
		//A[320 + i][(3 * 64) + ((64 + i - 0) % 64)] = 1; //x_3(0)
		//A[320 + i][(0 * 64) + ((64 + i - 46) % 64)] = 1; //d_1(46)
		//A[320 + i][(1 * 64) + ((64 + i - 45) % 64)] = 1; //e_0(45)
		//A[320 + i][(3 * 64) + ((64 + i - 44) % 64)] = 1; //e_2(44)
		//A[320 + i][(5 * 64) + ((64 + i - 45) % 64)] = 1; //e_5(45)

		//Seven equation
		//A[384 + i][(4 * 64) + ((64 + i - 0) % 64)] = 1; //x_4(0)
		//A[384 + i][(5 * 64) + ((64 + i - 0) % 64)] = 1; //x_5(0)
		//A[384 + i][(6 * 64) + ((64 + i - 0) % 64)] = 1; //x_6(0)
	}

	for (i = 0; i < 32; i++) {
		//Six equation (half) A1
		A[320 + i][(2 * 64) + ((64 + i + 55) % 64)] = 1; //x_2(55)
		A[320 + i][(3 * 64) + ((64 + i + 48) % 64)] = 1; //x_3(48)
		A[320 + i][(5 * 64) + ((64 + i + 23) % 64)] = 1; //x_5(23)
		A[320 + i][(1 * 64) + ((64 + i + 46) % 64)] = 1; //x_1(46)
		A[320 + i][(2 * 64) + ((64 + i + 54) % 64)] = 1; //x_2(54)
		A[320 + i][(4 * 64) + ((64 + i + 47) % 64)] = 1; //x_4(47)

		//Seven equation (half) A2
		//A[352 + i][(1 * 64) + ((64 + i + 3) % 64)] = 1; //x_1(3)
		//A[352 + i][(0 * 64) + ((64 + i + 20) % 64)] = 1; //x_0(20)
		//A[352 + i][(5 * 64) + ((64 + i + 41) % 64)] = 1; //x_5(41)
		//A[352 + i][(5 * 64) + ((64 + i + 43) % 64)] = 1; //x_5(43)
		//A[352 + i][(1 * 64) + ((64 + i + 2) % 64)] = 1; //x_1(2)
		//A[352 + i][(4 * 64) + ((64 + i + 3) % 64)] = 1; //x_4(3)

		//Seven equation (half) A3  +  21
		//A[352 + i][(1 * 64) + ((64 + i + 47 + 21) % 64)] = 1; //x_1(47)
		//A[352 + i][(2 * 64) + ((64 + i + 55 + 21) % 64)] = 1; //x_2(55)
		//A[352 + i][(4 * 64) + ((64 + i + 48 + 21) % 64)] = 1; //x_4(48)
		//A[352 + i][(0 * 64) + ((64 + i + 63 + 21) % 64)] = 1; //x_0(63)
		//A[352 + i][(1 * 64) + ((64 + i + 46 + 21) % 64)] = 1; //x_1(46)
		//A[352 + i][(5 * 64) + ((64 + i + 20 + 21) % 64)] = 1; //x_5(20)
		//A[352 + i][(5 * 64) + ((64 + i + 22 + 21) % 64)] = 1; //x_5(22)
	}

	

	/*****************************************************/
	//End of building linear equations
	/*****************************************************/

	GF2 det;
	DWORD Start, End;


#ifdef TEST_AVG
	Start = timeGetTime();
	for (int i = 0; i < 1000; i++)
	{
		solve(det, A, x, b); //solving the system of linear equation
	}
	End = timeGetTime();
	cout << "Cost of solving the system of linear equation: " << (double)(End - Start) / 1000 << "ms" << endl;
#else
	solve(det, A, x, b); //solving the system of linear equation
#endif // TEST_AVG


	vec_GF2 b_test, x_test;
	b_test.SetLength(2);
	x_test.SetLength(2);
	mat_GF2 A_test;
	clear(A_test);
	A_test.SetDims(2, 2);
	A_test[0][0] = A_test[1][1] = 1;
	A_test[0][1] = A_test[1][0] = 0;
	b[0] = 1;
	b[1] = 1;
	GF2 det_test;
	solve(det_test, A_test, x_test, b_test);
	cout << "det=" << det << endl;
	cout << "det_test=" << det_test << endl;
	cout << "x_test[0]=" << x_test[0] << " x_test[1]=" << x_test[1] << endl;
	//这里运行你的程序代码

	///*****************************************************/
	////Extracting the solution of the system of linear equation
	///*****************************************************/


	vector<vec_GF2> m_d;
	for (int i = 0; i < 17; i++)
	{
		m_d.push_back(vec_GF2());
		m_d[i].SetLength(64);
	}

	for (i = 0; i < 64; ++i) {
		m_d[0][i] = x[(0 * 64) + i];
		m_d[5][i] = x[(1 * 64) + i];
		m_d[10][i] = x[(2 * 64) + i];
		m_d[2][i] = x[(3 * 64) + i];
		m_d[7][i] = x[(4 * 64) + i];
		m_d[12][i] = x[(5 * 64) + i];
	}

	// 其他位置赋值

	for (int i = 0; i < 64; ++i) {
		m_d[6][i] = 1;
		m_d[11][i] = 1;
		m_d[15][i] = 1;
		m_d[16][i] = 1;
	}

	///*****************************************************/
	////End
	///*****************************************************/



#ifdef DEBUG
	cout << "The preimage is stored in the file named \"preimage\".\n";
#endif // DEBUG

	ofstream outfile;
	outfile.open(output_path, ios::out);
	if (!outfile) {
		cout << "Cannot open file \n";
		exit(1);
	}
#ifdef DEBUG
	debug_print_lane(m_d, 'd', 64, 17);
#endif // DEBUG

	output_print_lane(m_d, 'd', 64, 17, outfile);


	outfile.close();
}

vector<vec_GF2> get192zero(int bitsize, vector<vec_GF2>& hh)
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
	}

	//for padding

	e[11][0] = e[11][63] = e[16][0] = e[16][63] = 1;

	/*****************************************************/
	//End
	/*****************************************************/

	return e;
}



vector<vec_GF2> get_test_input(int bitsize, vector<vec_GF2>& hh)
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

	d[4][(64 - 27 + 43 + 64) % 64] = d[9][(64 - 27 + 43 + 64) % 64] = d[1][64 - 45] = d[16][64 - 45] = 1;

	/*****************************************************/
	//End
	/*****************************************************/

	return d;
}

vector<vec_GF2> get_test_input_2(int bitsize, vector<vec_GF2>& hh)
{
	int lanesize = bitsize / 64;
	vector<vec_GF2> e;
	for (int i = 0; i < 17; i++)
	{
		e.push_back(vec_GF2());
		e[i].SetLength(64);
	}


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