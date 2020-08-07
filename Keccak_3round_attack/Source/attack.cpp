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



void gaussian_init(int a[maxn][maxn], int x[maxn], int& equ, int& var)
{
	memset(a, 0, sizeof(a));
	memset(x, 0, sizeof(x));
	equ = 387, var = 384;
}


int gaussian_solve(int a[maxn][maxn], int x[maxn], int RC[5][maxn], int free_x[maxn], int& equ, int& var, int& free_num)
{
	int max_r, col, k;
	free_num = 0;
	for (k = 0, col = 0; k < equ && col < var; k++, col++) {
		max_r = k;
		for (int i = k + 1; i < equ; i++) {
			if (abs(a[i][col]) > abs(a[max_r][col]))
				max_r = i;
		}
		if (a[max_r][col] == 0) {
			k--;
			free_x[free_num++] = col;
			continue;
		}
		if (max_r != k) {
			for (int j = col; j < var + 1; j++)
				swap(a[k][j], a[max_r][j]);
		}
		for (int i = k + 1; i < equ; i++) {
			if (a[i][col] != 0) {
				for (int j = col; j < var + 1; j++)
					a[i][j] ^= a[k][j];
			}
		}
	}
	for (int i = k; i < equ; i++)
		if (a[i][col] != 0)
			return -1; //无解
	if (k < var)
		return var - k; //解不唯一，返回解个数
	//解唯一，生成解集
	for (int i = var - 1; i >= 0; i--) {
		x[i] = a[i][var];
		for (int j = i + 1; j < var; j++)
			x[i] ^= (a[i][j] && x[j]);
	}
	return 0;
}

int equ, var;
int a[maxn][maxn]; //增广矩阵
int RC[5][maxn];
int x[maxn]; //解集
int free_x[maxn];
int free_num;

void get_preimage(string input_path, string output_path)
{

	/*****************************************************/
	//Opening the file which contains the hash value.
	/*****************************************************/

	int i;

	vector<vec_GF2> h;
	const int hash_lane_num = HASH_BIT_SIZE / 64 + 1;
	rep(i, 0, hash_lane_num) {
		h.push_back(vec_GF2());
		h[i].SetLength(64);
	}


	// 做RC的预处理
	tKeccakLane rc[25];
	GetRoundConstants(rc);





	vector<vec_GF2> hh(h);

	iota_inverse(hh, 2, rc);
#ifdef DEBUG
	cout << "After Create RC and final hiota (32bits) ************\n";
	debug_print_lane(hh, 'i', 64, 5);
#endif // DEBUG






	/*****************************************************/
	//Building a system of linear equation. Ax = b where x is a vector containing lane variables.
	/*****************************************************/
	//// 抄一下RC获取值
	//GetRoundConstants(rc);
	// 转 int 数组
	rep(i, 0, 3)
	{
		cout << "RC[" << i << "]=";
		rep(j, 0, 64)
		{
			if (rc[i] & 1LL << j) {
				RC[i][j] = 1;
			}
			cout << RC[i][j];
		}
		cout << endl;
	}
	gaussian_init(a,x,equ,var);
	// b 赋值
	// for (i = 0; i < 64; ++i) {
	// 	// 注意求反当作常数
	// 	b[0 * 64 + i] = 0;
	// 	b[1 * 64 + i] = 0;
	// 	b[2 * 64 + i] = RC0[i];
	// 	b[3 * 64 + i] = 1;
	// 	b[4 * 64 + i] = 0;
	// }
	// for (i = 0; i < 32; ++i) {
	// 	b[320 + i] = RC1[i] + RC2[i] + 1;
	// 	//b[352 + i] = RC0[(64 + i + 20) % 64] + RC1[(64 + i + 20) % 64] + 1;
	// 	b[352 + i] = RC0[(64 + i + 63) % 64];
	// }

	// b赋值 （增广矩阵最右）
	rep(i, 0, 64)
	{
		a[0 * 64 + i][384] = 0;
		a[1 * 64 + i][384] = 0;
		a[2 * 64 + i][384] = RC[0][i];
		a[3 * 64 + i][384] = 1;
		a[4 * 64 + i][384] = 0;
	}
	rep(i, 0, 32)
	{
		if (i == 1) {
			a[320 + i][384] = RC[1][i] ^ RC[2][i] ^ 1;
		}
		else {
			a[320 + i][384] = RC[1][i] ^ RC[2][i] ^ 1;
		}
		// a[352 + i][384] = RC[0][(64 + i + 63 + 21) % 64];
	}
	rep(i, 0, 31)
	{
		// a[320 + i][384] = RC[1][i] ^ RC[2][i] ^ 1;
		if (i == 1) {
			a[352 + i][384] = RC[0][(64 + i + 63 + 21) % 64];
		}
		else {
			a[352 + i][384] = RC[0][(64 + i + 63 + 21) % 64];
		}
	}

	// A赋值
	rep(i, 0, 64)
	{
		//First
		a[0 + i][(0 * 64) + ((64 + i + 0) % 64)] = 1; //x_0(0)
		a[0 + i][(1 * 64) + ((64 + i + 0) % 64)] = 1; //x_1(0)
		a[0 + i][(2 * 64) + ((64 + i + 0) % 64)] = 1; //x_2(0)

		//Second
		a[64 + i][(3 * 64) + ((64 + i + 0) % 64)] = 1; //x_3(0)
		a[64 + i][(4 * 64) + ((64 + i + 0) % 64)] = 1; //x_4(0)
		a[64 + i][(5 * 64) + ((64 + i + 0) % 64)] = 1; //x_5(0)

		//Third
		a[128 + i][(0 * 64) + ((64 + i + 0) % 64)] = 1; //x_0(0)
		a[128 + i][(3 * 64) + ((64 + i + 2) % 64)] = 1; //x_3(2)
		a[128 + i][(5 * 64) + ((64 + i + 21) % 64)] = 1; //x_5(21)

		//Fourth
		a[192 + i][(1 * 64) + ((64 + i + 28) % 64)] = 1; //x_1(28)
		a[192 + i][(4 * 64) + ((64 + i + 58) % 64)] = 1; //x_4(58)

		//Fifth
		a[256 + i][(2 * 64) + ((64 + i + 61) % 64)] = 1; //x_2(61)
		a[256 + i][(5 * 64) + ((64 + i + 21) % 64)] = 1; //x_5(21)
	}
	rep(i, 0, 32)
	{
		//Sixth
		a[320 + i][(2 * 64) + ((64 + i + 55) % 64)] = 1; //x_2(55)
		a[320 + i][(3 * 64) + ((64 + i + 48) % 64)] = 1; //x_3(48)
		a[320 + i][(5 * 64) + ((64 + i + 23) % 64)] = 1; //x_5(23)
		a[320 + i][(1 * 64) + ((64 + i + 46) % 64)] = 1; //x_1(46)
		a[320 + i][(2 * 64) + ((64 + i + 54) % 64)] = 1; //x_2(54)
		a[320 + i][(4 * 64) + ((64 + i + 47) % 64)] = 1; //x_4(47)

		//Seventh 2rd
		// a[354 + i][(1 * 64) + ((64 + i + 47 + 21) % 64)] = 1; //x_1(47)
		// a[354 + i][(2 * 64) + ((64 + i + 55 + 21) % 64)] = 1; //x_2(55)
		// a[354 + i][(4 * 64) + ((64 + i + 48 + 21) % 64)] = 1; //x_4(48)
		// a[354 + i][(0 * 64) + ((64 + i + 63 + 21) % 64)] = 1; //x_0(63)
		// a[354 + i][(1 * 64) + ((64 + i + 46 + 21) % 64)] = 1; //x_1(46)
		// a[354 + i][(5 * 64) + ((64 + i + 20 + 21) % 64)] = 1; //x_5(20)
		// a[354 + i][(5 * 64) + ((64 + i + 22 + 21) % 64)] = 1; //x_5(22)
	}
	rep(i, 0, 31)
	{
		//Sixth
		// a[320 + i][(2 * 64) + ((64 + i + 55) % 64)] = 1; //x_2(55)
		// a[320 + i][(3 * 64) + ((64 + i + 48) % 64)] = 1; //x_3(48)
		// a[320 + i][(5 * 64) + ((64 + i + 23) % 64)] = 1; //x_5(23)
		// a[320 + i][(1 * 64) + ((64 + i + 46) % 64)] = 1; //x_1(46)
		// a[320 + i][(2 * 64) + ((64 + i + 54) % 64)] = 1; //x_2(54)
		// a[320 + i][(4 * 64) + ((64 + i + 47) % 64)] = 1; //x_4(47)

		//Serven 1st

		//Seventh 2nd
		a[352 + i][(1 * 64) + ((64 + i + 47 + 21) % 64)] = 1; //x_1(47)
		a[352 + i][(2 * 64) + ((64 + i + 55 + 21) % 64)] = 1; //x_2(55)
		a[352 + i][(4 * 64) + ((64 + i + 48 + 21) % 64)] = 1; //x_4(48)
		a[352 + i][(0 * 64) + ((64 + i + 63 + 21) % 64)] = 1; //x_0(63)
		a[352 + i][(1 * 64) + ((64 + i + 46 + 21) % 64)] = 1; //x_1(46)
		a[352 + i][(5 * 64) + ((64 + i + 20 + 21) % 64)] = 1; //x_5(20)
		a[352 + i][(5 * 64) + ((64 + i + 22 + 21) % 64)] = 1; //x_5(22)
	}
	a[383][0] = 1;
	a[384][66] = 1;

	/*****************************************************/
	//End of building linear equations
	/*****************************************************/

	GF2 det;
	DWORD Start, End;


#ifdef TEST_AVG
	Start = timeGetTime();
	for (int i = 0; i < 1000; i++)
	{
		//solve(det, A, x, b); //solving the system of linear equation
		int k = gaussian_solve();
}
	End = timeGetTime();
	cout << "Cost of solving the system of linear equation: " << (double)(End - Start) / 1000 << "ms" << endl;
#else
	//solve(det, A, x, b); //solving the system of linear equation
	int k = gaussian_solve(a, x, RC, free_x, equ, var, free_num);
#endif // TEST_AVG
	cout << "=== Answer ===" << endl;
	cout << "k=" << k << endl;
	rep(i, 0, 6)
	{
		cout << "x[" << i << "]=";
		rep(j, 0, 64)
		{
			cout << x[i * 64 + j];
		}
		cout << endl;
	}
	int temp;
	cin >> temp;

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