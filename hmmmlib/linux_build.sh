#!/bin/bash



# Run fresh cmake
fresh_cmake () {
    echo -e "Initializing fresh cmake build:"
    rm -r build
    mkdir build
    cd build
    cmake .. && echo -e "cmake: success\n"
}


# Run make
run_make () {
    make && echo -e "Built into build/\nmake: success\n"
}


# Test exe
test_exe () {
    echo -e "Executing executable:"
    ../build/exeHMMLIB && echo -e "executable: success\n"
}

# Test pinding
test_pinding () {
    echo -e "Running python binding test(s):"
    python ../../test_framework/test.py && echo -e "python binding tests: success"; 
}


fresh_cmake && run_make && test_exe && test_pinding