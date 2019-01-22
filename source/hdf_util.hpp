#ifndef HDF_UTIL
#define HDF_UTIL 1

#include <H5Cpp.h>
#include <string>
#include <sstream>
#include "hdsim.hpp"

using std::string;
using std::vector;
using std::pair;
using H5::Group;
using H5::PredType;
using H5::DataSpace;
using H5::DSetCreatPropList;
using H5::DataSet;
using H5::DataType;

/*! \brief Master function for writing vectors to hdf5 files
\param file Either an actual file or a group within a file
\param data Data to be written
\param caption Name of dataset
\param dt Data type
*/
template<class T> void write_std_vector_to_hdf5
(const Group& file,
	const vector<T>& data,
	const string& caption,
	const DataType& dt)
{
	hsize_t dimsf[1];
	dimsf[0] = static_cast<hsize_t>(data.size());
	DataSpace dataspace(1, dimsf);

	DSetCreatPropList plist;
	if (dimsf[0]>100000)
		dimsf[0] = 100000;
	plist.setChunk(1, dimsf);
	plist.setDeflate(6);

	DataSet dataset = file.createDataSet
		(H5std_string(caption),
			dt,
			dataspace,
			plist);
	dataset.write(&data[0], dt);
}

/*! \brief Writes floating point data to hdf5
\param file Either an actual file or a group within a file
\param data Data to be written
\param caption Name of dataset
*/
void write_std_vector_to_hdf5
(const Group& file,
	const vector<double>& data,
	const string& caption);

/*! \brief Writes integer data to hdf5
\param file Either an actual file or a group within a file
\param data Data to be written
\param caption Name of dataset
*/
void write_std_vector_to_hdf5
(const Group& file,
	const vector<int>& data,
	const string& caption);

/*! \brief Writes char data to hdf5
\param file Either an actual file or a group within a file
\param data Data to be written
\param caption Name of dataset
*/
void write_std_vector_to_hdf5
(const Group& file,
	const vector<unsigned char>& data,
	const string& caption);


//! \brief Container for snapshot data
class Snapshot
{
public:

	//! \brief Default constructor
	Snapshot(void);

	//! \brief Copy constructor
	//! \param source Source
	Snapshot(const Snapshot& source);

	//! \brief Mesh points
	vector<double> edges;

	//! \brief Computational cells
	vector<Primitive> cells;

	//! \brief Time
	double time;

	//! \brief Cycle number
	int cycle;

	//! \brief Total energy radiated away
	double Ecool;
};

/*! \brief Load snapshot data into memory
\param fname File name
\param mpioverride Flag for not reading mpi data when MPI is on
\return Snapshot data
*/
Snapshot read_hdf5_snapshot(const string& fname);


/*!
\brief Writes the simulation data into an HDF5 file
\param sim The hdsim class of the simulation
\param fname The name of the output file
\param appendices Additional data to be written to snapshot
\param a_names The name sof the appendices
*/
void write_snapshot_to_hdf5(hdsim & sim, string const& fname,std::vector<vector<double> > appendices = std::vector<vector<double> >(),
	std::vector<std::string> a_names = std::vector<std::string>());

string int2str(int n);

#endif // HDF_UTIL
