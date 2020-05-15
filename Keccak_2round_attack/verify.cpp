#include "verify.h"
#include "common.h"
using namespace std;

/*******  initialization functions ************/


void KeccakP1600_Initialize(void* state)
{
	memset(state, 0, 1600 / 8);
}

void KeccakP1600_StaticInitialize(void)
{
	if (sizeof(tKeccakLane) != 8) {
		printf("tKeccakLane should be 64-bit wide\n");

	}
	KeccakP1600_InitializeRoundConstants();
	KeccakP1600_InitializeRhoOffsets();
}


static int LFSR86540(UINT8* LFSR)
{
	int result = ((*LFSR) & 0x01) != 0;
	if (((*LFSR) & 0x80) != 0)
		/* Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1 */
		(*LFSR) = ((*LFSR) << 1) ^ 0x71;
	else
		(*LFSR) <<= 1;
	return result;
}

void KeccakP1600_InitializeRoundConstants(void)
{
	UINT8 LFSRstate = 0x01;
	unsigned int i, j, bitPosition;

	for (i = 0; i < maxNrRounds; i++) {
		KeccakRoundConstants[i] = 0;
		for (j = 0; j < 7; j++) {
			bitPosition = (1 << j) - 1; /* 2^j-1 */
			if (LFSR86540(&LFSRstate))
				KeccakRoundConstants[i] ^= (tKeccakLane)1 << bitPosition;
		}
	}
}

void KeccakP1600_InitializeRhoOffsets(void)
{
	unsigned int x, y, t, newX, newY;

	KeccakRhoOffsets[index(0, 0)] = 0;
	x = 1;
	y = 0;
	for (t = 0; t < 24; t++) {
		KeccakRhoOffsets[index(x, y)] = ((t + 1) * (t + 2) / 2) % 64;
		newX = (0 * x + 1 * y) % 5;
		newY = (2 * x + 3 * y) % 5;
		x = newX;
		y = newY;
	}
}




/******* 5 operation functions ************/


#define ROL64(a, offset) ((offset != 0) ? ((((tKeccakLane)a) << offset) ^ (((tKeccakLane)a) >> (64-offset))) : a)

static void theta(tKeccakLane* A)
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

static void rho(tKeccakLane* A)
{
	unsigned int x, y;

	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		A[index(x, y)] = ROL64(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

static void pi(tKeccakLane* A)
{
	unsigned int x, y;
	tKeccakLane tempA[25];

	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		tempA[index(x, y)] = A[index(x, y)];
	for (x = 0; x < 5; x++) for (y = 0; y < 5; y++)
		A[index(0 * x + 1 * y, 2 * x + 3 * y)] = tempA[index(x, y)];
}

static void chi(tKeccakLane* A)
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

static void iota(tKeccakLane* A, unsigned int indexRound)
{
	A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}




/*******  permutation functions ************/


void KeccakP1600Round(tKeccakLane* state, unsigned int indexRound)
{

	theta(state);
#ifdef DEBUG
	printf("$$$$$$$$$$$ After theta $$$$$$$$$$$$$\n");
	print1600state(state);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif // DEBUG

	rho(state);

#ifdef DEBUG
	printf("$$$$$$$$$$$ After rho $$$$$$$$$$$$$\n");
	print1600state(state);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif // DEBUG
	pi(state);

#ifdef DEBUG
	printf("$$$$$$$$$$$ After pi $$$$$$$$$$$$$\n");
	print1600state(state);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif // DEBUG
	chi(state);

#ifdef DEBUG
	printf("$$$$$$$$$$$ After chi $$$$$$$$$$$$$\n");
	print1600state(state);
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif // DEBUG

	iota(state, indexRound);

}

void KeccakF1600Permute(tKeccakLane* state, unsigned int nRounds) {
	for (int n = 0; n < nRounds; n++) {
		KeccakP1600Round(state, n);
	}
}



/****  Compute the hash value of input *****/

void Keccak(int reducedRound, unsigned int rate, unsigned int capacity, const UINT64* input,
	unsigned int inputBitLen, UINT64* output,
	unsigned int outputBitLen, int level) {

	if ((rate + capacity) != 1600 || (rate % 8) != 0) {
		printf("rate or capacity is illegal.\n");
		return;
	}


	tKeccakLane state[25];
	unsigned int blockSize = 0;
	unsigned int i;
	unsigned rateInWords = rate / 64;


	/* Initialize the state */
	KeccakP1600_StaticInitialize();
	KeccakP1600_Initialize(state);

	/* Absorb all the input blocks */
	int pos = 0;
	int blocknum = 0;

	while (inputBitLen > pos) {
		blocknum++;
		blockSize = MIN(inputBitLen - pos, rate);
		//state xor message block
		for (i = 0; i < blockSize; i++) {
			UINT64 mask = (UINT64)1 << ((pos + i) % 64);
			state[i / 64] ^= input[(pos + i) / 64] & mask;

		}
		pos += blockSize;
#ifdef DEBUG
		printf("************************ After XOR *******************************\n");
		print1600state(state);
		printf("******************************************************************\n\n");
#endif // DEBUG
		//no padding
		if (blockSize == rate) {
			KeccakF1600Permute(state, reducedRound);

#ifdef DEBUG
			printf("************** After Incircle Permutation (bsize=rate) *************\n");
			print1600state(state);
			printf("********************************************************************\n\n");
#endif // DEBUG
			blockSize = 0;
		}
	}

	/* Do the padding and switch to the squeezing phase */
	/* Add the first bit of padding */
	state[blockSize / 64] ^= (UINT64)1 << (blockSize % 64);


#ifdef DEBUG
	printf("************** After First Bit of Padding *************\n");
	print1600state(state);
	printf("*******************************************************\n\n");
#endif // DEBUG

	/* If the first bit of padding is at position rate - 1,
	we need a whole new block for the second bit of padding */
	if (blockSize == (rate - 1)) {
		KeccakF1600Permute(state, reducedRound);
	}

	/* Add the second bit of padding */
	state[rateInWords - 1] ^= (UINT64)1 << 63;

#ifdef DEBUG
	printf("************** After Scecond bit of Padding *************\n");
	print1600state(state);
	printf("*********************************************************\n\n");
#endif // DEBUG

	/* Switch to the squeezing phase */
	KeccakF1600Permute(state, reducedRound);


#ifdef DEBUG
	printf("************** After The Final Permutation *************\n");
	print1600state(state);
	printf("********************************************************\n\n");
#endif // DEBUG
	/* Squeeze out all the output blocks */
	pos = 0;
	while (outputBitLen - pos > 0) {
		blockSize = MIN(outputBitLen, rate);
		for (i = 0; i < blockSize; i++) {
			UINT64 mask = (UINT64)1 << (i % 64);
			output[(pos + i) / 64] ^= (state[i / 64] & mask) << (pos % 64);
		}
		pos += blockSize;
		if (outputBitLen - pos > 0) {
			KeccakF1600Permute(state, reducedRound);
			blockSize = 0;
		}
	}
}


/*** print the message and hash value to the screen ***/
void printHash(UINT64 message[], int inputLen, UINT64* hash, int outputLen) {
	printf("-The message (lenght = %d): \n", inputLen);
	for (int i = 0; i <= inputLen / 64; i++) {
		printf("%016llX ", message[i]);
	}
	printf("\n-The hash (length = %d): \n", outputLen);
	for (int i = 0; i <= (outputLen - 1) / 64; i++) {
		printf("%016llX ", hash[i]);
	}
	printf("\n");
}


/*** print the state ***/
void print1600state(tKeccakLane* state) {
	printf("==================STATE=====================\n");
	for (int i = 0; i < 86; i++) {
		printf("=");
	}
	printf("\n");
	for (int i = 4; i >= 0; --i) {
		printf("|");
		for (int j = 0; j < 5; ++j) {
			printf("%016llX|", state[i * 5 + j]);
		}
		printf("\n");
	}
	for (int i = 0; i < 86; i++) {
		printf("=");
	}
	printf("\n\n");
}

void score(UINT64* hash, int outputLen, int reducedRound) {
	const int table_delta[5] = { 1,2,3,10,20 };
	int ln = 0;

	for (int i = 0; i <= outputLen / 64; i++) {
		for (int j = 0; j < 64; j++) {
			if ((hash[i] >> j)& (1L)) {
				printf("-The length of 0's is %d for %d rounds Keccak256.\n\n", ln, reducedRound);
				if (ln < 64) {
					printf("*** length < 64 : Score = 0. ***\n");
				}
				else {
					printf("*** Score_%d = %d * %d = %d. ***\n", reducedRound, ln, table_delta[reducedRound - 1], ln * table_delta[reducedRound - 1]);
				}

				return;
			}
			else {
				ln++;
			}
		}
	}
}

UINT64 stringtoull(string str_64) {
	// 2进制字符串 转UINT64
	return std::strtoull(str_64.c_str(), NULL, 2);
}


void test_preimage(string input_path) {
#define		MESSAGE_BITLEN (1088*2-64)
	//#define  MESSAGE_BITLEN (1088-64)

	UINT64 message_block[MESSAGE_BITLEN / 64] = { 0 };
	unsigned int inputBitLen1 = MESSAGE_BITLEN;
	int reducedRound = 1;

	ifstream infile;
	infile.open(input_path);
	rep(i, 0, MESSAGE_BITLEN / 64) {
		string se;
		//	sd = vecGF2tostring(d[i],64);
		getline(infile, se);
		cout << "Readline [" << i << "]: " << se << endl;
		//	cout << "sd=" << sd << endl << "se=" << se << endl;
		//	// 是否需要反转
		//	reverse(sd.begin(), sd.end());
		reverse(se.begin(), se.end());
		message_block[i] = stringtoull(se);
	}

	cout << "***************PRINT MESSAGE*************\n";
	rep(i, 0, MESSAGE_BITLEN / 64) {
		printf("i=%d  %016llX\n", i, message_block[i]);
	}

	unsigned int outputBitLen = 256;
	UINT64* hash_1 = (UINT64*)calloc(outputBitLen / 64 + 1, sizeof(UINT64));

	Keccak(reducedRound, 1088, 512, message_block, inputBitLen1, hash_1, outputBitLen, 0);

	printf("M_1:\n");
	printHash(message_block, inputBitLen1, hash_1, outputBitLen);

	//compute score
	score(hash_1, outputBitLen, reducedRound);

}

//
//int main(int argc, const char* argv[]) {
//
//	// example1
//	UINT64 message_text[7] = { 0x00,0x00,0x00,0x00,0x00,0x01,0x00 };
//	unsigned int inputBitLen1 = 447;
//	int reducedRound = 1;
//
//	// example2
//	UINT64 message_1[5] = { 0x01,0x8000000000000000,0x00,0x00,0x00 };
//	unsigned int inputBitLen1 = 319;
//	int reducedRound = 1;
//
//	unsigned int outputBitLen = 256;
//	UINT64* hash_1 = (UINT64*)calloc(outputBitLen / 64 + 1, sizeof(UINT64));
//
//	Keccak(reducedRound, 1088, 512, message_1, inputBitLen1, hash_1, outputBitLen, 0);
//
//	printf("M_1:\n");
//	printHash(message_1, inputBitLen1, hash_1, outputBitLen);
//
//	//compute score
//	score(hash_1, outputBitLen, reducedRound);
//
//	return 0;
//}




