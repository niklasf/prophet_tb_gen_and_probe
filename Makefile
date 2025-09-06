
#-NDEBUG
# -fno-exceptions
# -O3
fast:
	clang++ -g -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64  -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full   -c -o main.o main.cpp
	clang++ -g -o main bitboard.o main.o uci.o -m64 -lpthread  -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64  -DNDEBUG -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full


build:
	clang++ -g -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64   -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full   -c -o bitboard.o bitboard.cpp
	clang++ -g -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64  -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full   -c -o uci.o uci.cpp
	clang++ -g -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64  -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full   -c -o main.o main.cpp
	clang++ -g -o main bitboard.o main.o uci.o -m64 -lpthread  -Wall -Wcast-qual -fno-exceptions -std=c++17  -pedantic -Wextra -Wshadow -Wmissing-declarations -m64  -DNDEBUG -O2 -funroll-loops -DIS_64BIT -DUSE_POPCNT     -flto=full

