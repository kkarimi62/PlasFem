module load gcc/4.8.2_gcc-4.4.6
c++ -fopenmp -I/home/kkarimi/opt/eigen-eigen-07105f7124f9 -I/home/kkarimi/opt/boost_1_59_0 -O2 -std=c++11 -static-libstdc++ *.cpp -o fem_run -w
