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
perf_link_flags := -O2 $(common_link_flags)
cov_link_flags := -O0 --coverage $(common_link_flags) $(common_test_link_flags)

main_object := main
prog_object_names := geometry kd_tree model_loader raycaster material image\
                     input renderer
test_common_object_names := test_utils
test_unit_object_names := geometry_test kd_tree_test utils_test \
                          model_loader_test raycaster_test lighting_utils_test \
                          image_test material_test vector_test input_test
test_perf_object_names := performance_main performance_utils \
                          sort_performance kd_performance

bin_objects = $(patsubst %, bin/%.o, $(main_object) $(prog_object_names))

test_prog_objects := $(patsubst %, testing/bin/%.o, $(prog_object_names))
test_common_objects := $(patsubst %, testing/bin/%.o, $(test_common_object_names))
test_unit_objects := $(patsubst %, testing/bin/%.o, $(test_unit_object_names))
test_perf_objects := $(patsubst %, testing/bin/%.o, $(test_perf_object_names))

cov_prog_objects = $(patsubst %, testing/cov/%.o, $(prog_object_names))
cov_common_objects = $(patsubst %, testing/cov/%.o, $(test_common_object_names))
cov_unit_objects = $(patsubst %, testing/cov/%.o, $(test_unit_object_names))

all_unit_objects := $(test_prog_objects) $(test_common_objects) $(test_unit_objects)
all_perf_objects := $(test_prog_objects) $(test_common_objects) $(test_perf_objects)
all_cov_objects := $(cov_prog_objects) $(cov_common_objects) $(cov_unit_objects)

all: bin/main testing/bin/unit testing/bin/performance testing/cov/main

bin/main : $(bin_objects)
	$(CC) $(bin_link_flags) -o bin/main $(bin_objects)

testing/bin/unit : $(all_unit_objects)
	$(CC) $(test_link_flags) -o testing/bin/unit $(all_unit_objects)

testing/bin/performance : $(all_perf_objects)
	$(CC) $(perf_link_flags) -o testing/bin/performance $(all_perf_objects)

testing/cov/main: $(all_cov_objects)
	$(CC) $(cov_link_flags) -o testing/cov/main $(all_cov_objects)

# for $$(notdir %) in the rule prerequisites
.SECONDEXPANSION:

bin/%.o: %.cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(bin_compile_flags) -c -MMD $< -o $@

$(test_prog_objects): %.o: $$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

$(test_common_objects): %.o: testing/$$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

$(test_unit_objects): %.o: testing/unit/$$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

$(test_perf_objects): %.o: testing/performance/$$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(test_compile_flags) -c -MMD $< -o $@

$(cov_prog_objects): %.o: $$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c -MMD $< -o $@

$(cov_common_objects): %.o: testing/$$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c -MMD $< -o $@ 

$(cov_unit_objects): %.o: testing/unit/$$(notdir %).cpp
	@mkdir -p $(@D)
	$(CC) $(common_compile_flags) $(cov_compile_flags) -c -MMD $< -o $@ 


-include $(wildcard bin/*.d)
-include $(wildcard testing/bin/*.d)
-include $(wildcard testing/cov/*.d)

.PHONY: clean
clean :
	rm -rf ./testing/cov
	rm -rf ./testing/bin
	rm -rf ./bin

.PHONY: test
test :
	./testing/bin/unit

.PHONY: coverage
coverage :
	@mkdir -p testing/cov
	-./testing/cov/main
	gcovr -r . --html --html-details -o testing/cov/coverage.html

.PHONY: performance
performance :
	./testing/bin/performance
