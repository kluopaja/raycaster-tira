CC := g++
common_compile_flags := -std=c++17 -Wall -Weffc++ -Wextra -Wsign-conversion
test_compile_flags := -O2
cov_compile_flags := -O0 --coverage

common_test_link_flags = -lgtest -lgtest_main -lgmock
test_link_flags := -O2 $(common_test_link_flags)
cov_link_flags := -O0 --coverage $(common_test_link_flags)

objects := geometry.o kd_tree.o
test_objects := geometry_test.o kd_tree_test.o utils_test.o test_utils.o

test_depends := $(patsubst %.o, test/%.d, $(test_objects) $(objects))
# cov_depends := $(patsubst %.cpp, cov/%.d, $(test_objects) $(objects))


all: cov/main test/main

test/main : $(patsubst %, test/%, $(test_objects) $(objects))
	@mkdir -p test
	$(CC) $(test_link_flags) -o test/main $(cov_objects)

cov/main : $(patsubst %, cov/%, $(test_objects) $(objects))
	@mkdir -p cov
	$(CC) $(cov_link_flags) -o cov/main $(cov_objects)

test/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

cov/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c -MMD $< -o $@ 



# -include $(test_depends)
# -include $(cov_depends)

.PHONY: clean
clean :
	rm -rf ./cov
	rm -rf ./test
	rm -rf ./bin

.PHONY: test
test :
	test/main

.PHONY: coverage
coverage :
	@mkdir -p cov
	cov/main
	gcovr -r . --html --html-details -o cov/coverage.html
