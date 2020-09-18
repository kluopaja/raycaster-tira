CC := g++
common_compile_flags := -std=c++17 -Wall -Weffc++ -Wextra -Wsign-conversion
cov_compile_flags := -O2 --coverage
cov_link_flags := -lgtest -lgtest_main -lgmock --coverage
objects := bin/geometry.o bin/kd_tree.o bin/utils.o
cov_objects := cov/geometry_test.o cov/geometry.o cov/test_utils.o \
               cov/kd_tree.o cov/kd_tree_test.o cov/utils_test.o


cov/test : $(cov_objects)
	@mkdir -p cov
	$(CC) $(cov_link_flags) -o cov/test $(cov_objects)

cov/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c $< -o $@ 

clean :
	rm -rf ./cov

coverage :
	@mkdir -p cov
	cov/test
	gcovr -r . --html --html-details -o cov/coverage.html
