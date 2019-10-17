SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
SOURCES := $(RAW_SOURCES)
LIB_FILE = libLag.a
CC := g++
CCC := gcc
ARCHIVER_FUNC := ar
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

set_environ_vars.sh: | external_libraries/include/H5Cpp.h external_libraries/boost_dump/boost_1_66_0/boost/container/static_vector.hpp
	$(eval MY_BOOST_PATH=`pwd`/external_libraries/boost_dump/boost_1_66_0)
	$(eval MY_HDF5_PATH=`pwd`/external_libraries/include)
	echo export\ CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH):$(MY_BOOST_PATH):$(MY_HDF5_PATH) > set_environ_vars.sh
	echo export\ HDF5_LIB_PATH=`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ LD_PATH=$(LD_PATH):`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ RICH_ROOT=`pwd` >> set_environ_vars.sh

external_libraries/include/H5Cpp.h: external_libraries/hdf5_dump/hdf5-1.10.5/c++/src/H5Cpp.h
	cd external_libraries/hdf5_dump/hdf5-1.10.5 && \
	./configure --enable-cxx --prefix=`cd ../.. && pwd`
	cd external_libraries/hdf5_dump/hdf5-1.10.5 && make
	cd external_libraries/hdf5_dump/hdf5-1.10.5 && make check
	cd external_libraries/hdf5_dump/hdf5-1.10.5 && make install

external_libraries/hdf5_dump/hdf5-1.10.5/c++/src/H5Cpp.h: | external_libraries/hdf5_dump/hdf5-1.10.5.tar.gz
	cd external_libraries/hdf5_dump/ && tar xf ./hdf5-1.10.5.tar.gz

external_libraries/boost_dump/boost_1_66_0/boost/container/static_vector.hpp: | external_libraries/boost_dump/boost_1_66_0.tar.bz2
	cd external_libraries/boost_dump/ && tar xf ./boost_1_66_0.tar.bz2

external_libraries/hdf5_dump/hdf5-1.10.5.tar.gz:
	mkdir -p external_libraries/hdf5_dump
	cd external_libraries/hdf5_dump && \
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz

external_libraries/boost_dump/boost_1_66_0.tar.bz2:
	mkdir -p external_libraries/boost_dump
	cd external_libraries/boost_dump && \
	wget 'http://sourceforge.net/projects/boost/files/boost/1.66.0/boost_1_66_0.tar.bz2/download'
	cd external_libraries/boost_dump && mv download boost_1_66_0.tar.bz2
