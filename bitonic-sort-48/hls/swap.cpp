#include "swap.hpp"


void swap1(PFOutputObj &data1, PFOutputObj &data2) {
	//if (data1.range(BIT_LOW,BIT_HIGH) < data2.range(BIT_LOW,BIT_HIGH)) {
	if (data1.hwPt < data2.hwPt) {
		std::swap(data1, data2);
	}
}

void swap2(PFOutputObj &data1, PFOutputObj &data2) {
	//if (data1.range(BIT_LOW,BIT_HIGH) > data2.range(BIT_LOW,BIT_HIGH)) {
	if (data1.hwPt > data2.hwPt) {
		std::swap(data1, data2);
	}
}
