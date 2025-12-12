
# CC = /opt/homebrew/bin/g++-15
CC = g++
flags =  -Wall -Wextra  -std=c++17 -fopenmp -O2
repair:
	$(CC) -g -c  block_compress_egtb.cpp repair.cpp encoder.cpp huffman.cpp $(flags)
	$(CC) -g -o block_compress_egtb.out block_compress_egtb.o repair.o encoder.o huffman.o $(flags)