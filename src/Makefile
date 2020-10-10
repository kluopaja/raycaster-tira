enable_parallel = 0

CC := g++
common_compile_flags := -std=c++17 -Wall -Wextra -Wsign-conversion
bin_compile_flags := -O2 -D ENABLE_PARALLEL=$(enable_parallel)
test_compile_flags := -O2
cov_compile_flags := -O0 --coverage

common_link_flags =  -lassimp -pthread
common_test_link_flags = -lgtest -lgtest_main -lgmock

bin_link_flags := -O2 $(common_link_flags)
ifeq ($(enable_parallel), 1)
  bin_link_flags += -ltbb
endif
test_link_flags := -O2 $(common_link_flags) $(common_test_link_flags)
cov_link_flags := -O0 --coverage $(common_link_flags) $(common_test_link_flags)

main_object := main
common_object_names := geometry kd_tree model_loader raycaster material image
test_object_names := geometry_test kd_tree_test utils_test test_utils \
                     model_loader_test raycaster_test lighting_utils_test \
                     image_test material_test


bin_objects = $(patsubst %, bin/%.o, $(main_object) $(common_object_names))
test_objects = $(patsubst %, test/%.o, $(test_object_names) \
               $(common_object_names))
cov_objects = $(patsubst %, cov/%.o, $(test_object_names) \
              $(common_object_names))

bin_depends := $(patsubst %, bin/%.d, $(main_object) $(common_object_names))
test_depends := $(patsubst %, test/%.d, $(test_object_names) \
                $(common_object_names))
cov_depends := $(patsubst %, cov/%.d, $(test_objects_names) \
               $(common_objects_names))

all: bin/main test/main cov/main

bin/main : $(bin_objects)
	@mkdir -p bin
	$(CC) $(bin_link_flags) -o bin/main $(bin_objects)

test/main : $(test_objects)
	@mkdir -p test
	$(CC) $(test_link_flags) -o test/main $(test_objects)

cov/main : $(cov_objects)
	@mkdir -p cov
	$(CC) $(cov_link_flags) -o cov/main $(cov_objects)

bin/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(bin_compile_flags) -c -MMD $< -o $@

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
	-cov/main
	gcovr -r . --html --html-details -o cov/coverage.html
