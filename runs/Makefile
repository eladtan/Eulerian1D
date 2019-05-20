CDIR = `pwd`
SDIR = $(shell dirname $(CDIR))/source
CC := g++
DUMP_FOLDER := temp_$(TARGET)
PREXECUTE :=
HDF5_LIB_PATH := /home/elad/huji-rich/external_libraries/lib
BOOST_DIR := /home/elad/huji-rich/external_libraries/boost_dump/boost_1_66_0
HDF_DIR := /home/elad/huji-rich/external_libraries/include
OPTIMIZATION_FLAGS := -I $(SDIR) -I $(BOOST_DIR) -I $(HDF_DIR)
ENDEMIC_SOURCES := $(shell find $(TARGET) -name '*.cpp')
ENDEMIC_OBJECTS := $(patsubst $(TARGET)/%.cpp,$(DUMP_FOLDER)/%.o,$(ENDEMIC_SOURCES))
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS += -O0 -g -pg -std=c++11
	LFLAGS := -pg
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS += -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS :=
	LFLAGS :=
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel)
	CCC := icc
	CC := mpiicpc
	OPTIMIZATION_FLAGS += -DRICH_MPI -O3 -ipo -xHost -fp-model precise -std=c++11 -DOMPI_SKIP_MPICXX
	LINT_FLAGS = 
	ARCHIVER_FUNC := xiar
else
	MODE := production
	OPTIMIZATION_FLAGS += -O3 -std=c++11
	LFLAGS :=
endif
CFLAGS := $(LINT_FLAGS) $(OPTIMIZATION_FLAGS)

$(DUMP_FOLDER)/normal_termination.res: $(DUMP_FOLDER)/test.exe

 $(DUMP_FOLDER)/test.exe: $(ENDEMIC_OBJECTS) | ../library_$(MODE)/libLag.a 	
	echo "$(ENDEMIC_SOURCES)"
	$(CC) $^ $| -o $@ $(LFLAGS) -lhdf5 -lhdf5_cpp -L $(HDF5_LIB_PATH)

../library_$(MODE)/librich.a:
	$(MAKE) -C .. MODE=$(MODE)

$(ENDEMIC_OBJECTS): $(DUMP_FOLDER)/%.o: $(TARGET)/%.cpp
	export LD_LIBRARY_PATH=$(HDF5_LIB_PATH):$(LD_LIBRARY_PATH)
	mkdir -p $(DUMP_FOLDER)
	@echo $(LD_LIBRARY_PATH)
	$(CC) $(CFLAGS) -c $< -o $@ -I ..
