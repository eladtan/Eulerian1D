#ifndef MPI_COMM_HPP
#define MPI_COMM_HPP 1
#ifdef RICH_MPI
#include <mpi.h>
#include "Extensive.hpp"
#include "Primitive.hpp"
#include <vector>
#include <array>
#include "ExactRS.hpp"

#define NGHOSTCELLS 2

std::array<Primitive, NGHOSTCELLS * 2> SendRecvPrimitive(std::vector<Primitive> const& cells);

std::array<double, NGHOSTCELLS * 2> SendRecvEdges(std::vector<double> const& edges);

void RedistributeExtensives(std::vector<Extensive> &cells, std::vector<double> &edges, std::vector<Primitive> &pcells);

void ConsolidateData(std::vector<Primitive> &cells, std::vector<double> &edges, std::vector<std::vector<
	double> > & append, double &Ecool);

#endif

#endif //MPI_COMM_HPP