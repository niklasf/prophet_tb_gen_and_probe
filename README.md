# Prophet TB: Generator and Probing

This repository contains code to generate and probe chess distance-to-mate (DTM) endgame tablebases.

I was able to generate the tablebase for up to 6 pieces with a Ryzen 9 5950x CPU (16 cores 32 hyper-threads) and 96GB RAM in 14 days.

The 6-men tablebase is block-compressed with 32768 entries per block and is 1111 GiB large.

For now, probing requires mmaping the compressed files and is thus only available on unix-based systems.

A single DTM query takes <0.1ms on my system, where the files are stored on a M.2 SSD.

In the full tablebase, there are two files for each piece configuration, e.g. in KQKNP the side-to-move has a queen against a knight and pawn and in KNPKQ the side-to-move has knight and pawn against queen. (There is just one file for symmetric configurations, e.g. KRPKRP).

To save space you may choose to only store one table per piece configuration, which brings size down to 449 GiB.
However, since all moves have to be expanded and probed for positions where the table is now missing, probing takes on average 10 times longer.

## Building

First, you need to build [zstd](https://github.com/facebook/zstd).

```
git submodule update --init
cd zstd
make
cd ..
```


Run
```
make prophet
```
to build the shard library `build/prophet.so` used for probing.

For a demo, build
```
make mates
```
and run `build/longest_mates.out [threads] [lines]` which finds the longest forced mate for each piece configuration with `[threads]` and by substituting `[lines]` with `1` it also finds the forced move sequence (otherwise set `[lines] = 0`). 


For generating the tablebase, build
```
make gen
```
and run `build/gen_main.out [max_npiece] [threads]` to compute with `[threads]` number of threads the tables for all piece configuration with a maximum number of `[max_npiece]` pieces.
Results are stored in `egtbs/`

On my 32 thread system, 4-men takes <1 minute and 5-men takes <2 hours to compute. (5-men needs 7 GiB disk space and 16GB RAM).