#include "hdf_util.hpp"
#include <cassert>
#ifdef RICH_MPI
#include "mpi_comm.hpp"
#endif

using namespace H5;

void write_std_vector_to_hdf5
(const Group& file,
	const vector<double>& data,
	const string& caption)
{
	FloatType datatype(PredType::NATIVE_DOUBLE);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
		(file,
			data,
			caption,
			datatype);
}

void write_std_vector_to_hdf5
(const Group& file,
	const vector<unsigned char>& data,
	const string& caption)
{
	PredType datatype(PredType::NATIVE_UCHAR);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
	(file,
		data,
		caption,
		datatype);
}


void write_std_vector_to_hdf5
(const Group& file,
	const vector<int>& data,
	const string& caption)
{
	IntType datatype(PredType::NATIVE_INT);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
		(file,
			data,
			caption,
			datatype);
}


Snapshot::Snapshot(void) :
	edges(),
	cells(),
	time(),
	cycle(),
	Ecool(){}

Snapshot::Snapshot(const Snapshot& source) :
	edges(source.edges),
	cells(source.cells),
	time(source.time),
	cycle(source.cycle),
	Ecool(source.Ecool){}

namespace 
{
	template<class T> vector<T> read_vector_from_hdf5
		(const Group& file,
			const string& caption,
			const DataType& datatype)
	{
		DataSet dataset = file.openDataSet(caption);
		DataSpace filespace = dataset.getSpace();
		hsize_t dims_out[2];
		filespace.getSimpleExtentDims(dims_out, NULL);
		const size_t NX = static_cast<size_t>(dims_out[0]);
		vector<T> result(NX);
		dataset.read(&result[0], datatype);
		return result;
	}

	vector<double> read_double_vector_from_hdf5
		(Group& file, string const& caption)
	{
		return read_vector_from_hdf5<double>
			(file,
				caption,
				PredType::NATIVE_DOUBLE);
	}

	vector<int> read_int_vector_from_hdf5
		(const Group& file,
			const string& caption)
	{
		return read_vector_from_hdf5<int>
			(file,
				caption,
				PredType::NATIVE_INT);
	}

	vector<unsigned char> read_char_vector_from_hdf5
	(const Group& file,
		const string& caption)
	{
		return read_vector_from_hdf5<unsigned char>
			(file,
				caption,
				PredType::NATIVE_UCHAR);
	}
}


void write_snapshot_to_hdf5(hdsim &sim, string const& fname, std::vector<vector<double> > appendices,
	std::vector<std::string> a_names)
{
#ifdef RICH_MPI
	// consolidate data
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<Primitive> oldcells = sim.GetCells();
	std::vector<double> oldedges = sim.GetEdges();
	std::vector<vector<double> > oldappend = appendices;
	double OldEcool = sim.GetEcool();
	ConsolidateData(sim.GetCells(), sim.GetEdges(),appendices,sim.GetEcool());
	if (rank == 0)
	{
#endif
		H5File file(H5std_string(fname), H5F_ACC_TRUNC);
		Group geometry = file.createGroup("/geometry");
		Group hydrodynamic = file.createGroup("/hydrodynamic");

		// General
		write_std_vector_to_hdf5
		(file,
			vector<double>(1, sim.GetTime()),
			"time");
		write_std_vector_to_hdf5
		(file,
			vector<int>(1, sim.GetCycle()),
			"cycle");
		write_std_vector_to_hdf5
		(file,
			vector<double>(1, sim.GetEcool()),
			"Ecool");

		// Geometry  
		write_std_vector_to_hdf5
		(geometry, sim.GetEdges(), "edges");

		// Hydrodynamic
		vector<Primitive> const& cells = sim.GetCells();
		size_t N = cells.size();
		vector<double> density(N), pressure(N), velocity(N), LastCool(N);
		std::vector<unsigned char> stickers(N);
		for (size_t i = 0; i < N; ++i)
		{
			density[i] = cells[i].density;
			pressure[i] = cells[i].pressure;
			velocity[i] = cells[i].velocity;
			stickers[i] = cells[i].sticker;
			LastCool[i] = cells[i].LastCool;
		}

		write_std_vector_to_hdf5
		(hydrodynamic, density,
			"density");
		write_std_vector_to_hdf5
		(hydrodynamic,
			pressure,
			"pressure");
		write_std_vector_to_hdf5
		(hydrodynamic,
			velocity,
			"velocity");
		write_std_vector_to_hdf5
		(hydrodynamic,
			LastCool,
			"LastCool");
		write_std_vector_to_hdf5
		(hydrodynamic,
			stickers,
			"stickers");

		// appendices
		if (a_names.size() > 0)
		{
			Group Gappend = file.createGroup("/appendices");
			assert(a_names.size() == appendices.size());
			size_t Nappend = appendices.size();
			for (size_t i = 0; i < Nappend; ++i)
				write_std_vector_to_hdf5(Gappend, appendices[i], a_names[i]);
		}
#ifdef RICH_MPI
		sim.GetCells() = oldcells;
		sim.GetEdges() = oldedges;
		appendices = oldappend;
		sim.GetEcool() = OldEcool;
		sim.GetCells().shrink_to_fit();
		sim.GetEdges().shrink_to_fit();
		for (size_t i = 0; i < appendices.size(); ++i)
			appendices[i].shrink_to_fit();
	}
#endif
}

Snapshot read_hdf5_snapshot
(const string &fname)
{
	Snapshot res;
	H5File file(fname, H5F_ACC_RDONLY);
	Group g_geometry = file.openGroup("geometry");
	Group g_hydrodynamic = file.openGroup("hydrodynamic");

	// Mesh points
	res.edges = read_double_vector_from_hdf5(g_geometry, "edges");

	// Hydrodynamic
	{
		const vector<double> density =
			read_double_vector_from_hdf5(g_hydrodynamic, "density");
		const vector<double> pressure =
			read_double_vector_from_hdf5(g_hydrodynamic, "pressure");
		const vector<double> velocity =
			read_double_vector_from_hdf5(g_hydrodynamic, "velocity");
		const vector<unsigned char> stickers =
			read_char_vector_from_hdf5(g_hydrodynamic, "stickers");
		res.cells.resize(density.size());
		for (size_t i = 0; i<res.cells.size(); ++i) 
		{
			res.cells.at(i).density = density.at(i);
			res.cells.at(i).pressure = pressure.at(i);
			res.cells.at(i).velocity = velocity.at(i);
			res.cells.at(i).sticker = stickers.at(i);
		}
	}

	// Misc
	{
		const vector<double> time =
			read_double_vector_from_hdf5(file, "time");
		res.time = time.at(0);
		const vector<double> Ecool =
			read_double_vector_from_hdf5(file, "Ecool");
		res.Ecool = Ecool.at(0);
		const vector<int> cycle =
			read_int_vector_from_hdf5(file, "cycle");
		res.cycle = cycle.at(0);
	}
	return res;
}

string int2str(int n)
{
	stringstream ss;
	ss << n;
	return ss.str();
}