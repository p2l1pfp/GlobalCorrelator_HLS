#!/usr/bin/env python3
# -*- coding: utf-8 -*-def gen():

import os
import sys

DATA_SIZE = 64

def gen(filename, pipeline_II):

	with open(filename, 'a') as f:

		main_s = ''
		f.write('#include "swap.hpp"\n\n')

		block = 2
		while (block <= DATA_SIZE):
			step = int(block / 2)
			while (step >= 1):

				s = 'void block%d_step%d_net(PFOutputObj datas[DATA_SIZE]) {\n' % (block, step)
				f.write(s)

				main_s = main_s + '	block%d_step%d_net(datas);\n' % (block, step)

				for idx1 in range(0, DATA_SIZE):
					idx2 = idx1 ^ step
					if (idx1 >= idx2):
						continue
					if ((idx1 & block) != 0):
						s = '	swap1(datas[%d], datas[%d]);\n' % (idx1, idx2)
						f.write(s)
					else:
						s = '	swap2(datas[%d], datas[%d]);\n' % (idx1, idx2)
						f.write(s)

				s = '}\n\n'
				f.write(s)

				step = int(step / 2)
			block = block * 2

		f.write('\nvoid sorting_network(PFOutputObj datas[DATA_SIZE]) {\n')
		f.write('#pragma HLS pipeline II=%d rewind\n\n' % (pipeline_II))
		f.write(main_s)
		f.write('}\n')


def gen_header(filename_header):
	with open(filename_header, 'a') as f:
		f.write('#ifndef __SORTNING_NETWORK_HPP__\n')
		f.write('#define __SORTNING_NETWORK_HPP__\n\n')
		f.write('#include "bitonic_sort.hpp"\n\n')
		f.write('void sorting_network(PFOutputObj datas[DATA_SIZE]);\n')
		f.write('\n#endif\n')


if __name__ == '__main__':

	argvs = sys.argv
	argc = len(argvs)

	if (argc != 2):
		print('Usage: # python3 %s pipeline_II' % argvs[0])
		quit()

	pipeline_II = int(argvs[1])

	filename = './sorting_network.cpp'
	filename_header = './sorting_network.hpp'
	print('Generate Sorting Network')

	gen(filename, pipeline_II)
	gen_header(filename_header)
