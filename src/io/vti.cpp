#include "vti.hpp"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace io {

static std::string header_vti(int i0, int i1, int j0, int j1, const VTIOptions &opts) {
	std::ostringstream oss;
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <ImageData WholeExtent=\"" << i0 << " " << i1 << " "
			<< j0 << " " << j1 << " 0 0\" Origin=\"" 
			<< opts.origin[0] << " " << opts.origin[1] << " " << opts.origin[2]
			<< "\" Spacing=\"" << opts.spacing[0] << " " << opts.spacing[1] << " "
			<< opts.spacing[2] << "\">\n";
	oss << "    <Piece Extent=\"" << i0 << " " << i1 << " "
			<< j0 << " " << j1 << " 0 0\">\n";
	if (opts.pointData) {
		oss << "      <PointData Scalars=\"" << opts.arrayName << "\">\n";
		oss << "        <DataArray type=\"Float64\" Name=\"" << opts.arrayName
				<< "\" format=\"ascii\"/>\n";
		oss << "      </PointData>\n";
		oss << "      <CellData/>\n";
	} else {
		oss << "      <PointData/>\n";
		oss << "      <CellData Scalars=\"" << opts.arrayName << "\">\n";
		oss << "        <DataArray type=\"Float64\" Name=\"" << opts.arrayName
				<< "\" format=\"ascii\"/>\n";
		oss << "      </CellData>\n";
	}
	oss << "    </Piece>\n";
	oss << "  </ImageData>\n";
	oss << "</VTKFile>\n";
	return oss.str();
}

static std::string header_pvti(int gi0, int gi1, int gj0, int gj1, const VTIOptions &opts, int nranks) {
	std::ostringstream oss;
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <PImageData WholeExtent=\"" << gi0 << " " << gi1 << " "
			<< gj0 << " " << gj1 << " 0 0\" Origin=\"" 
			<< opts.origin[0] << " " << opts.origin[1] << " " << opts.origin[2]
			<< "\" Spacing=\"" << opts.spacing[0] << " " << opts.spacing[1] << " "
			<< opts.spacing[2] << "\">\n";
	if (opts.pointData) {
		oss << "    <PPointData Scalars=\"" << opts.arrayName << "\">\n";
		oss << "      <PDataArray type=\"Float64\" Name=\"" << opts.arrayName << "\"/>\n";
		oss << "    </PPointData>\n";
		oss << "    <PCellData/>\n";
	} else {
		oss << "    <PPointData/>\n";
		oss << "    <PCellData Scalars=\"" << opts.arrayName << "\">\n";
		oss << "      <PDataArray type=\"Float64\" Name=\"" << opts.arrayName << "\"/>\n";
		oss << "    </PCellData>\n";
	}
	// pieces added by caller
	return oss.str();
}

void write_vti_parallel_2d(const std::string &stem,
													 MPI_Comm comm,
													 int global_nx, int global_ny,
													 int local_x0, int local_y0,
													 int local_nx, int local_ny,
													 const double *local_data,
													 const VTIOptions &opts) {
	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	// Local extents in point indexing
	const int i0 = local_x0;
	const int i1 = local_x0 + local_nx - 1;
	const int j0 = local_y0;
	const int j1 = local_y0 + local_ny - 1;

	// Overlapped extents with ghost layers (clamped to domain)
	const int g = std::max(0, opts.ghostLevel);
	const int li0 = std::max(0, i0 - g);
	const int li1 = std::min(global_nx - 1, i1 + g);
	const int lj0 = std::max(0, j0 - g);
	const int lj1 = std::min(global_ny - 1, j1 + g);

	// Write per-rank .vti
		{
			std::ostringstream fn;
			fn << stem << "_" << rank << ".vti";
			std::ofstream out(fn.str());
			out << "<?xml version=\"1.0\"?>\n";
			out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
			out << "  <ImageData WholeExtent=\"" << li0 << " " << li1 << " " << lj0 << " " << lj1 << " 0 0\" ";
		out << "Origin=\"" << opts.origin[0] << " " << opts.origin[1] << " " << opts.origin[2] << "\" ";
		out << "Spacing=\"" << opts.spacing[0] << " " << opts.spacing[1] << " " << opts.spacing[2] << "\">\n";
			out << "    <Piece Extent=\"" << li0 << " " << li1 << " " << lj0 << " " << lj1 << " 0 0\">\n";
		if (opts.pointData) {
			out << "      <PointData Scalars=\"" << opts.arrayName << "\">\n";
			out << "        <DataArray type=\"Float64\" Name=\"" << opts.arrayName << "\" format=\"ascii\" NumberOfComponents=\"1\" >\n";
		} else {
			out << "      <PointData/>\n";
			out << "      <CellData Scalars=\"" << opts.arrayName << "\">\n";
			out << "        <DataArray type=\"Float64\" Name=\"" << opts.arrayName << "\" format=\"ascii\" NumberOfComponents=\"1\" >\n";
		}
		out << std::setprecision(17);
			for (int gj = lj0; gj <= lj1; ++gj) {
				for (int gi = li0; gi <= li1; ++gi) {
					// Clamp to this piece's interior to synthesize ghost values by edge replication
					const int ci = std::clamp(gi, i0, i1) - i0;
					const int cj = std::clamp(gj, j0, j1) - j0;
					const double v = local_data[cj * local_nx + ci];
				out << v << " ";
			}
			out << "\n";
		}
		out << "        </DataArray>\n";
		if (opts.pointData) {
			out << "      </PointData>\n";
			out << "      <CellData/>\n";
		} else {
			out << "      </CellData>\n";
		}
		out << "    </Piece>\n";
		out << "  </ImageData>\n";
		out << "</VTKFile>\n";
		out.close();
	}

	// Gather piece metadata to rank 0 and write .pvti
		struct Piece { int i0,i1,j0,j1,rank; } local{ li0,li1,lj0,lj1,rank };
	std::vector<Piece> pieces;
	if (rank == 0) pieces.resize(size);
	MPI_Gather(&local, sizeof(Piece), MPI_BYTE,
						 rank==0? pieces.data():nullptr, sizeof(Piece), MPI_BYTE,
						 0, comm);

	if (rank == 0) {
		std::ostringstream fn;
		fn << stem << ".pvti";
		std::ofstream out(fn.str());
		out << "<?xml version=\"1.0\"?>\n";
			out << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
			out << "  <PImageData WholeExtent=\"0 " << (global_nx-1) << " 0 " << (global_ny-1) << " 0 0\" ";
		out << "Origin=\"" << opts.origin[0] << " " << opts.origin[1] << " " << opts.origin[2] << "\" ";
			out << "Spacing=\"" << opts.spacing[0] << " " << opts.spacing[1] << " " << opts.spacing[2] << "\" ";
			out << "GhostLevel=\"" << std::max(0, opts.ghostLevel) << "\">\n";
		if (opts.pointData) {
			out << "    <PPointData Scalars=\"" << opts.arrayName << "\">\n";
			out << "      <PDataArray type=\"Float64\" Name=\"" << opts.arrayName << "\" NumberOfComponents=\"1\"/>\n";
			out << "    </PPointData>\n";
			out << "    <PCellData/>\n";
		} else {
			out << "    <PPointData/>\n";
			out << "    <PCellData Scalars=\"" << opts.arrayName << "\">\n";
			out << "      <PDataArray type=\"Float64\" Name=\"" << opts.arrayName << "\" NumberOfComponents=\"1\"/>\n";
			out << "    </PCellData>\n";
		}
		for (const auto &p : pieces) {
			out << "    <Piece Extent=\"" << p.i0 << " " << p.i1 << " " << p.j0 << " " << p.j1 << " 0 0\" ";
			out << "Source=\"" << stem << "_" << p.rank << ".vti\"/>\n";
		}
		out << "  </PImageData>\n";
		out << "</VTKFile>\n";
		out.close();
	}
}

} // namespace io

