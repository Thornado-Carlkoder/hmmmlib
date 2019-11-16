#!/bin/bash



# Run fresh cmake
rm -r build
mkdir build
cd build
cmake .. && echo -e "cmake successfull\n"



# Run make
make && echo -e "Built into build/\nmake successfull"

# Test exe
echo -e "\nExecuting executable:"
../build/exeHMMLIB

# Test pinding
echo -e "\nRunning python binding test(s):"
python ../../test_framework/test.py; 
