
# CC = /opt/homebrew/bin/g++-15
CC = g++
flags =  -Wall -Wextra  -std=c++17 -fopenmp -O3
ldeflate = -I libdeflate -L libdeflate/build -l deflate
deflate:
	$(CC) -g -o compress.out block_compress_egtb.cpp $(flags) $(ldeflate)