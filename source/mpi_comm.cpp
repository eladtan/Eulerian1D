#include "mpi_comm.hpp"
#include <array>
#include <assert.h>
#include <iostream>
#include <cmath>
#ifdef RICH_MPI

namespace
{
	std::array<double, 7> Primitive2Array(Primitive const& cell)
	{
		std::array<double, 7> res;
		res[0] = cell.density;
		res[1] = cell.energy;
		res[2] = cell.entropy;
		res[3] = cell.LastCool;
		res[4] = cell.velocity;
		res[5] = cell.pressure;
		res[6] = static_cast<double>(cell.sticker);
		return res;
	}

	Primitive Array2Primitive(std::array<double, 7> const& data)
	{
		Primitive res;
		res.density = data[0];
		res.energy = data[1];
		res.entropy = data[2];
		res.LastCool = data[3];
		res.velocity = data[4];
		res.pressure = data[5];
		res.sticker = static_cast<unsigned char>(data[6]);
		return res;
	}

	std::array<double, 5> Extensive2Array(Extensive const& cell)
	{
		std::array<double, 5> res;
		res[0] = cell.mass;
		res[1] = cell.energy;
		res[2] = cell.entropy;
		res[3] = cell.momentum;
		res[4] = cell.et;
		return res;
	}

	Extensive Array2Extensive(std::array<double, 5> const& data)
	{
		Extensive res;
		res.mass = data[0];
		res.energy = data[1];
		res.entropy = data[2];
		res.momentum = data[3];
		res.et = data[4];
		return res;
	}

	std::vector<double> ExtensiveVec2Double(std::vector<Extensive> const& cells, std::vector<Primitive> const& pcells)
	{
		size_t N = cells.size();
		std::vector<double> res;
		res.reserve(7 * N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 5> temp = Extensive2Array(cells[i]);
			for (size_t j = 0; j < 5; ++j)
				res.push_back(temp[j]);
			res.push_back(pcells[i].LastCool);
			res.push_back(static_cast<double>(pcells[i].sticker));
		}
		return res;
	}

	std::vector<Extensive> VecDouble2Extensive(std::vector<double> const& data, std::vector<Primitive> &cells)
	{
		size_t N = data.size() / 7;
		if (N * 7 != data.size())
			assert(false);
		std::vector<Extensive> res(N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 5> temp;
			for (size_t j = 0; j < 5; ++j)
				temp[j] = data[i * 7 + j];
			res[i] = Array2Extensive(temp);
			cells[i].LastCool = data[i * 7 + 5];
			cells[i].sticker = static_cast<unsigned char>(data[i * 7 + 6]);
		}
		return res;
	}

	std::vector<double> PrimitiveVec2Double(std::vector<Primitive> const& cells)
	{
		size_t N = cells.size();
		std::vector<double> res;
		res.reserve(7 * N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 7> temp = Primitive2Array(cells[i]);
			for (size_t j = 0; j < 7; ++j)
				res.push_back(temp[j]);
		}
		return res;
	}

	std::vector<Primitive> VecDouble2Primitive(std::vector<double> const& data)
	{
		size_t N = data.size() / 7;
		assert(data.size() % 7 == 0);
		std::vector<Primitive> res(N);
		for (size_t i = 0; i < N; ++i)
		{
			std::array<double, 7> temp;
			for (size_t j = 0; j < 7; ++j)
				temp[j] = data[i * 7 + j];
			res[i] = Array2Primitive(temp);
		}
		return res;
	}
}

std::array<Primitive, NGHOSTCELLS * 2> SendRecvPrimitive(std::vector<Primitive> const& cells)
{
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::array<double, 7 * NGHOSTCELLS> data, data2, recv;
	std::array<double, 7> temp;
	std::array<Primitive, NGHOSTCELLS * 2> res;
	MPI_Request req, req2;
	assert(cells.size() > 2);
	if (rank == 0)
	{
		size_t N = cells.size();
		for (size_t i = 0; i < NGHOSTCELLS; ++i)
		{
			temp = Primitive2Array(cells[N - NGHOSTCELLS + i]);
			std::copy(temp.begin(), temp.end(), data.begin() + i * 7);
		}
		MPI_Isend(&data[0], 7 * NGHOSTCELLS, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &req);
		MPI_Recv(&recv[0], 7 * NGHOSTCELLS, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (size_t i = 0; i < NGHOSTCELLS; ++i)
		{
			std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
			res[NGHOSTCELLS + i] = Array2Primitive(temp);
		}
		MPI_Wait(&req, MPI_STATUS_IGNORE);
	}
	else
	{
		if (rank == (ws - 1))
		{
			for (size_t i = 0; i < NGHOSTCELLS; ++i)
			{
				temp = Primitive2Array(cells[i]);
				std::copy(temp.begin(), temp.end(), data.begin() + i * 7);
			}
			MPI_Isend(&data[0], 7 * NGHOSTCELLS, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
			MPI_Recv(&recv[0], 7 * NGHOSTCELLS, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (size_t i = 0; i < NGHOSTCELLS; ++i)
			{
				std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
				res[i] = Array2Primitive(temp);
			}
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Status status;
			size_t N = cells.size();
			for (size_t i = 0; i < NGHOSTCELLS; ++i)
			{
				temp = Primitive2Array(cells[N - NGHOSTCELLS + i]);
				std::copy(temp.begin(), temp.end(), data2.begin() + i * 7);
			}
			for (size_t i = 0; i < NGHOSTCELLS; ++i)
			{
				temp = Primitive2Array(cells[i]);
				std::copy(temp.begin(), temp.end(), data.begin() + i * 7);
			}
			MPI_Isend(&data[0], 7 * NGHOSTCELLS, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
			MPI_Isend(&data2[0], 7 * NGHOSTCELLS, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req2);
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 7 * NGHOSTCELLS, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
				{
					std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
					res[i] = Array2Primitive(temp);
				}
			}
			else
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
				{
					std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
					res[NGHOSTCELLS + i] = Array2Primitive(temp);
				}
			}
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], 7 * NGHOSTCELLS, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
				{
					std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
					res[i] = Array2Primitive(temp);
				}
			}
			else
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
				{
					std::copy(recv.begin() + 7 * i, recv.begin() + 7 * (i + 1), temp.begin());
					res[NGHOSTCELLS + i] = Array2Primitive(temp);
				}
			}
			MPI_Wait(&req, MPI_STATUS_IGNORE);
			MPI_Wait(&req2, MPI_STATUS_IGNORE);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

std::array<double, NGHOSTCELLS * 2> SendRecvEdges(std::vector<double> const & edges)
{
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::array<double, NGHOSTCELLS * 2> res;
	std::array<double, NGHOSTCELLS> recv;
	MPI_Request req, req2;
	assert(edges.size() > 3);
	if (rank == 0)
	{
		size_t N = edges.size() - 1;
		MPI_Isend(&edges[N - NGHOSTCELLS], NGHOSTCELLS, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &req);
		MPI_Recv(&recv[0], NGHOSTCELLS, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (size_t i = 0; i < NGHOSTCELLS; ++i)
			res[NGHOSTCELLS + i] = recv[i];
		MPI_Wait(&req, MPI_STATUS_IGNORE);
	}
	else
	{
		if (rank == (ws - 1))
		{
			MPI_Isend(&edges[1], NGHOSTCELLS, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req);
			MPI_Recv(&recv[0], NGHOSTCELLS, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (size_t i = 0; i < NGHOSTCELLS; ++i)
				res[i] = recv[i];
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Status status;
			size_t N = edges.size() - 1;
			MPI_Isend(&edges[1], NGHOSTCELLS, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &req);
			MPI_Isend(&edges[N - NGHOSTCELLS], NGHOSTCELLS, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req2);
			MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], NGHOSTCELLS, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
					res[i] = recv[i];
			}
			else
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
					res[NGHOSTCELLS + i] = recv[i];
			}
			MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&recv[0], NGHOSTCELLS, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_SOURCE == (rank - 1))
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
					res[i] = recv[i];
			}
			else
			{
				for (size_t i = 0; i < NGHOSTCELLS; ++i)
					res[NGHOSTCELLS + i] = recv[i];
			}
			MPI_Wait(&req, MPI_STATUS_IGNORE);
			MPI_Wait(&req2, MPI_STATUS_IGNORE);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

void RedistributeExtensives(std::vector<Extensive> &cells, std::vector<double> &edges, std::vector<Primitive> &pcells)
{
	int nlocal = cells.size();
	int ntotal = 0;
	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	// Do we need ro rebalance?
	double maxload = (1.0*nlocal * ws) / ntotal;
	int newload = 0;
	double res_factor = std::pow(static_cast<double>(ntotal) / (ws * 1000), -0.33333);
	res_factor *= 0.05;
	res_factor = std::min(res_factor, 0.15);
	if (maxload > res_factor || maxload < (1-res_factor))
		newload = 1;
	int shouldcalc = 0;
	MPI_Allreduce(&newload, &shouldcalc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (shouldcalc == 0)
		return;
	if (rank == 0)
		std::cout << "Redistribute data" << std::endl;
	size_t index_lower = (rank * ntotal) / ws;
	size_t index_upper = ((rank + 1) * ntotal) / ws - 1;
	std::vector<int> nperproc(ws, 0), disp(ws, 0), newn(ws, 0);
	nlocal *= 7;
	assert(nlocal > 0);
	MPI_Gather(&nlocal, 1, MPI_INT, &nperproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<double> tosend = ExtensiveVec2Double(cells, pcells);
	assert(tosend.size() == nlocal);
	std::vector<double> torecv(4);
	if (rank == 0)
	{
		torecv.resize(ntotal * 7);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(&tosend[0], tosend.size(), MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0], MPI_DOUBLE, 0,
		MPI_COMM_WORLD);
	tosend.resize(7 * (index_upper - index_lower + 1));
	// send the data
	if (rank == 0)
	{
		for (size_t i = 0; i < ws; ++i)
			newn[i] = 7 * (((i + 1) * ntotal) / ws - (i * ntotal) / ws);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = newn[i - 1] + disp[i - 1];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatterv(&torecv[0], &newn[0], &disp[0], MPI_DOUBLE, &tosend[0],
		7 * (index_upper - index_lower + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	pcells.resize(tosend.size() / 7);
	cells = VecDouble2Extensive(tosend, pcells);
	// deal with the edges
	if (rank == 0)
	{
		torecv.resize((ntotal + 1));
		for (size_t i = 0; i < ws; ++i)
		{
			nperproc[i] /= 7;
		}
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	size_t Nedges = edges.size() - 1;
	tosend.resize(Nedges);
	for (size_t i = 0; i < Nedges; ++i)
	{
		tosend[i] = edges[i + 1];
	}
	MPI_Gatherv(&tosend[0], Nedges, MPI_DOUBLE, &torecv[1], &nperproc[0], &disp[0],
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		torecv[0] = edges[0];
	}
	tosend.resize((index_upper - index_lower + 2));
	if (rank == 0)
	{
		for (size_t i = 0; i < ws; ++i)
		{
			newn[i] /= 7;
			newn[i]++;
			if (i > 0)
				disp[i] = newn[i - 1] + disp[i - 1] - 1;
		}
	}
	size_t Nnew = (index_upper - index_lower + 2);
	MPI_Scatterv(&torecv[0], &newn[0], &disp[0], MPI_DOUBLE, &tosend[0], Nnew, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	edges.resize(Nnew);
	for (size_t i = 0; i < Nnew; ++i)
	{
		edges[i] = tosend[i];
	}
	// clear memory
	edges.shrink_to_fit();
	cells.shrink_to_fit();
	pcells.shrink_to_fit();
}

void ConsolidateData(std::vector<Primitive>& cells, std::vector<double>& edges, std::vector<std::vector<double> >& append, double &Ecool)
{
	int nlocal = cells.size();
	int ntotal = 0;
	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int rank = 0, ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	std::vector<int> nperproc(ws, 0), disp(ws, 0);
	nlocal *= 7;
	MPI_Gather(&nlocal, 1, MPI_INT, &nperproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::vector<double> tosend = PrimitiveVec2Double(cells);
	std::vector<double> torecv(2);
	assert(nlocal > 0);
	if (rank == 0)
	{
		torecv.resize(ntotal * 7);
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&tosend[0], tosend.size(), MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0], MPI_DOUBLE, 0,
		MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		cells = VecDouble2Primitive(torecv);
	// deal with the edges
	if (rank == 0)
	{
		torecv.resize(ntotal + 1);
		for (size_t i = 0; i < nperproc.size(); ++i)
			nperproc[i] /= 7;
		for (size_t i = 1; i < ws; ++i)
			disp[i] = nperproc[i - 1] + disp[i - 1];
	}
	MPI_Gatherv(&edges[1], edges.size() - 1, MPI_DOUBLE, &torecv[1], &nperproc[0], &disp[0],
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		torecv[0] = edges[0];
		edges = torecv;
	}
	// deal with appendices
	for (size_t i = 0; i < append.size(); ++i)
	{
		nlocal = append[i].size();
		ntotal = 0;
		MPI_Gather(&nlocal, 1, MPI_INT, &nperproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank == 0)
		{
			ntotal = nlocal;
			torecv.resize(ntotal);
			disp[0] = 0;
			for (size_t i = 1; i < ws; ++i)
			{
				disp[i] = nperproc[i - 1] + disp[i - 1];
				ntotal += nperproc[i];
			}
			torecv.resize(ntotal);
		}
		MPI_Gatherv(&append[i][0], nlocal, MPI_DOUBLE, &torecv[0], &nperproc[0], &disp[0],
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0)
			append[i] = torecv;
	}
	double NewTot = 0;
	MPI_Reduce(&Ecool, &NewTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
		Ecool = NewTot;
	MPI_Barrier(MPI_COMM_WORLD);
}
#endif