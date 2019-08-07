SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
SOURCES := $(RAW_SOURCES)
LIB_FILE = libLag.a
CC := g++
CCC := gcc
ARCHIVER_FUNC := ar
BOOST_DIR := /home/elad/huji-rich/external_libraries/boost_dump/boost_1_66_0
HDF_DIR := /home/elad/huji-rich/external_libraries/include
HDF5_LIB_PATH := /home/elad/huji-rich/external_libraries/lib
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg -std=c++11
	LINT_FLAGS :=
else ifeq ($(MODE),intel)
	CC := icpc
	CCC := icc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS :=
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel)
	CCC := icc
	CC := mpiicpc
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -ipo -xHost -fp-model precise -std=c++11 -DOMPI_SKIP_MPICXX
	LINT_FLAGS = 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel_habanero)
	CCC := icc
	CC := mpiicpc
	OPTIMIZATION_FLAGS += -DRICH_MPI -O3 -ipo -xHost -fp-model precise -std=c++11 -DOMPI_SKIP_MPICXX
	LINT_FLAGS = 
	ARCHIVER_FUNC := xiar
	HDF5_LIB_PATH := /rigel/spice/users/es3640/Rich/external_libraries/lib
	BOOST_DIR := /rigel/spice/users/es3640/Rich/external_libraries/boost_dump/boost_1_66_0
	HDF_DIR := /rigel/spice/users/es3640/Rich/external_libraries/include
else
	MODE = production
	OPTIMIZATION_FLAGS := -O3 -march=native -std=c++11 -fno-expensive-optimizations
endif
OPTIMIZATION_FLAGS += -I $(HDF_DIR) -I $(BOOST_DIR) -I $(SOURCE_DIR)
LIBRARY_FOLDER := library_$(MODE)
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(SOURCES))

$(LIBRARY_FOLDER)/$(LIB_FILE): $(OBJECTS)
	$(ARCHIVER_FUNC) cr $@ $^

$(OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CC) -c $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $@
	$(CC) -MM $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$@,g' -i $(LIBRARY_FOLDER)/$*.d

-include $(OBJECTS:.o=.d)

clean:
	rm -rf ./$(LIBRARY_FOLDER)

set_environ_vars.sh: 
	echo export\ LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(HDF5_LIB_PATH) >> set_environ_vars.sh
