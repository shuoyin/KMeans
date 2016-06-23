#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/polar.h>
#include <data/xmipp_fftw.h>
#include <data/histogram.h>
#include <data/numerical_tools.h>
#include <data/xmipp_program.h>
#include <data/filters.h>
#include <reconstruction/fourier_filter.h>
#include <vector>
#include <queue>
#include <string>
#include <time.h>

#include "mpi.h"
#include "NetSim.h"

double fitTwoImages(MultidimArray<double> &fixed,MultidimArray<double> &I, double sigma, Matrix2D<double> &M, bool reverse=false);

double fit(MultidimArray<double> &fixed, MultidimArray<double> &I, double sigma, Matrix2D<double> &M);

void readImage(MetaData &md, Image<double> &I, size_t objId, bool applyGeo);

MPI::Comm &comm=MPI::COMM_WORLD;
int id, numprocs;

class Program{
public:
	int K0, KN;
	int Nimgs;
	FileName fin;
	MetaData md;
	size_t Xdim, Ydim;
	double sigma;
	MultidimArray<double> corrM;
	MultidimArray<double> adjancent;
public:
	Program(std::string s="", int k0=3, int kn=4);
	bool run();
};

class Cluster{
public:
	MultidimArray<double> thisP;
	MultidimArray<double> nextP;
	MultidimArray<double> simM;
	MultidimArray<double> netSimM;
	std::vector<int> belongings;
	std::vector<int> newbelongings;
	std::vector<double> corrV;
public:
	Cluster(){};
	void update(MetaData &md, double sigma);
	double getNetSim(int index);
};

class KMeans{
public:
	int K;
	int iter;
	int Nimgs;
	FileName name;
	MetaData md;
	int Xdim, Ydim;
	double sigma;
	std::vector<Cluster> centers;
	std::vector<int> oldclusters, newclusters;
public:
	KMeans(int k, int it, FileName fn, int n, int xdim, int yidm, double s):K(k), iter(it), name(fn), Nimgs(n), Xdim(xdim), Ydim(yidm), sigma(s){md.read(name);}
	void initialize(MultidimArray<double> &corrM);
	void createCluster(Cluster &c, MultidimArray<double> &corrM);
	void splitCluster(Cluster &node, Cluster &node1, Cluster &node2, MultidimArray<double> &corrM);
	bool run(MultidimArray<double> &corrM);
	void runDivisive(int finalN,MultidimArray<double> &corrM);
	void writeResult();
};

void printCenterInfo(KMeans &km, std::string name);