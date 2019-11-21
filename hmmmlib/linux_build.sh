#!/bin/bash



# Run fresh cmake
fresh_cmake () {
    rm -r build
    mkdir build
    cd build
    cmake .. && echo -e "cmake successfull\n"
}


# Run make
run_make () {
    make && echo -e "Built into build/\nmake successful"
}


# Test exe
test_exe () {
    echo -e "\nExecuting executable:"
    ../build/exeHMMLIB
}

# Test pinding
test_pinding () {
    echo -e "\nRunning python binding test(s):"
    python ../../test_framework/test.py; 
}


fresh_cmake && run_make && test_exe && test_pinding