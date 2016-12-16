Windows and *nix tested.
Compiler: C++11
For *nix systems, use gcc-4.7
g++ -std=c+11 main.cpp -o hamper

Four parameters are expected:
1. the full path and file name of the target genome
2. the full path and file name of the reads file
3. word size (note, a default maximum of 10 is set. this is the expected maximum for 16GB RAM machines. This can be altered within main.cpp)
4. threshold between 0.0 and 1.0 (0%-100%) - value of % ID required for a read to be considered a match
