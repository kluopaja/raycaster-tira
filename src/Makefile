CC := g++
common_compile_flags := -std=c++17 -Wall -Wextra -Wsign-conversion
test_compile_flags := -O2
cov_compile_flags := -O0 --coverage

common_link_flags =  -lassimp
common_test_link_flags = -lgtest -lgtest_main -lgmock
test_link_flags := -O2 $(common_link_flags) $(common_test_link_flags)
cov_link_flags := -O0 --coverage $(common_link_flags) $(common_test_link_flags)

object_names := geometry kd_tree model_loader raycaster
test_object_names := geometry_test kd_tree_test utils_test test_utils \
                     model_loader_test raycaster_test

test_objects = $(patsubst %, test/%.o, $(test_object_names) $(object_names))
cov_objects = $(patsubst %, cov/%.o, $(test_object_names) $(object_names))

test_depends := $(patsubst %, test/%.d, $(test_object_names) $(object_names))
cov_depends := $(patsubst %, cov/%.d, $(test_objects_names) $(objects_names))

all: cov/main test/main

test/main : $(test_objects)
	@mkdir -p test
	$(CC) $(test_link_flags) -o test/main $(test_objects)

cov/main : $(cov_objects)
	@mkdir -p cov
	$(CC) $(cov_link_flags) -o cov/main $(cov_objects)

test/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

cov/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c -MMD $< -o $@ 

-include $(test_depends)
-include $(cov_depends)

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
