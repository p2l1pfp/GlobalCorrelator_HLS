#include "bitonic_sort.hpp"

#include <cstdio>
#include <algorithm>


// Number of test case
#define TEST_NUM 32

struct TestData {
	PFOutputObj datas[DATA_SIZE];
	PFOutputObj datas_acutual[DATA_SIZE];
	PFOutputObj datas_expected[DATA_SIZE];
};


void makeTestData(TestData& td);
int check(TestData& td);
int tb_bitonic_sort();

bool range_sort(PFOutputObj data1, PFOutputObj data2) {
	//if (data1.range(BIT_LOW,BIT_HIGH) > data2.range(BIT_LOW,BIT_HIGH)) {
	if (data1.hwPt > data2.hwPt) {
            return true;
	} else {
            return false;
        }
}

int main() {
	// return number of error
	return tb_bitonic_sort();
}

typedef PFOutputObj OBJ;

int tb_bitonic_sort() {

	hls::stream<OBJ> axis_in[DATA_SIZE];
	hls::stream<OBJ> axis_out[DATA_SIZE];

	// Make test data
	TestData testdatas[TEST_NUM];
	for (int id = 0; id < TEST_NUM; ++id) {
		makeTestData(testdatas[id]);
	}

	// Set input data
	for (int id = 0; id < TEST_NUM; ++id) {
		for (int i = 0; i < DATA_SIZE; ++i) {
			axis_in[i].write(testdatas[id].datas[i]);
		}
	}

	// Execute
	for (int id = 0; id < TEST_NUM; ++id) {
		bitonic_sort(axis_in, axis_out);
	}

	// Get output
	bool flag = true;
	for (int id = 0; id < TEST_NUM; ++id) {
		for (int i = 0; i < DATA_SIZE; ++i) {
			testdatas[id].datas_acutual[i] = axis_out[i].read();
		}
	}

	// Check
	int err = 0;
	for (int id = 0; id < TEST_NUM; ++id) {
		err += check(testdatas[id]);
	}
	return err;

}


void makeTestData(TestData& td) {
	const int MAX = 2047;
	const int MIN = -2047;
	for (int i = 0; i < DATA_SIZE; ++i) {
            if (i < DATA_NONZERO) {
		td.datas[i].hwPt = rand()% (MAX - MIN + 1) + MIN;
		td.datas_expected[i] = td.datas[i];
            } else {
		td.datas[i].hwPt = 0;
		td.datas_expected[i] = td.datas[i];
            }
	}
	std::sort(td.datas_expected, td.datas_expected + DATA_NONZERO, range_sort);
}


int check(TestData& td) {
	int ret = 0;
	for (int i = 0; i < DATA_SIZE; ++i) {
                std::cout<<i<<"\t: actual = "<<td.datas_acutual[i].hwPt<<"\t\texpected = "<<td.datas_expected[i].hwPt<<std::endl;
		if (td.datas_acutual[i].hwPt != td.datas_expected[i].hwPt)
			ret++;
	}
	return ret;
}
