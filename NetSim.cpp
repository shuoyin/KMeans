#include "NetSim.h"

using namespace std;

void CreateAdjacentMatrix(MultidimArray<double> &SimIn, MultidimArray<double> &Out, int k, const Cmp cmp)
{
	Out.initZeros(SimIn);
	for(int i=0; i<SimIn.ydim; i++){
		priority_queue<Elem, vector<Elem>, Cmp> q(cmp);
		for(int j=0; j<SimIn.xdim; j++){
			if(j==i) continue;
			q.push(Elem(j, dAij(SimIn, i, j)));
		}
		for(int j=0; j<k; j++){
			dAij(Out, i, q.top().index)=1;
			// dAij(Out, q.top().index, i)=1;
			q.pop();
		}
	}
	MultidimArray<double> tmp=Out;
	for(int i=0; i<SimIn.ydim; i++){
		for(int j=0; j<SimIn.xdim; j++){
			dAij(Out, i, j) = dAij(tmp, i, j)*dAij(tmp, j, i);
		}
	}
}

void CreateAdjacentMatrix(vector<double> &SimIn, vector<double> &Out, int k, const Cmp cmp)
{
	int n=SimIn.size();
	Out.resize(n,0);
	priority_queue<Elem, vector<Elem>, Cmp> q(cmp);
	for(int i=0; i<n; i++)
		q.push(Elem(i, SimIn[i]));
	for(int i=0; i<k; i++){
		Out[q.top().index]=1;
		q.pop();
	}
}

void SNN(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> sharedNear(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, sharedNear, k, cmp);

	for(int i=0; i<SimOut.ydim; i++){
		for(int j=0; j<SimOut.xdim; j++){
			if(i==j){
				dAij(SimOut, i, j) = k;
				continue;
			}
			if(j<i){
				dAij(SimOut, i, j) = dAij(SimOut, j, i);
				continue;
			}
			dAij(SimOut, i, j) = 0;
			for(int m=0; m<sharedNear.xdim; m++)
				dAij(SimOut, i, j)+=dAij(sharedNear, i, m)*dAij(sharedNear, j, m);
		}
	}
}

void JaccardIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	for(int i=0; i<SimOut.ydim; i++){
		for(int j=0; j<SimOut.xdim; j++){
			if(i==j){
				dAij(SimOut, i, j) = 1;
				continue;
			}
			if(j<i){
				dAij(SimOut, i, j) = dAij(SimOut, j, i);
				continue;
			}
			double intersection=0;
			double Union=0;
			for(int m=0; m<adjacent.xdim; m++){
				intersection+=dAij(adjacent, i, m)*dAij(adjacent, j, m);
				Union+=max(dAij(adjacent, i, m), dAij(adjacent, j, m));
			}
			// Union-=intersection;
			dAij(SimOut, i, j) = intersection/Union;
		}
	}
}

void SorensenIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	for(int i=0; i<SimOut.ydim; i++){
		for(int j=0; j<SimOut.xdim; j++){
			if(i==j){
				dAij(SimOut, i, j) = 1;
				continue;
			}
			if(j<i){
				dAij(SimOut, i, j) = dAij(SimOut, j, i);
				continue;
			}
			double num=0;
			double den=0;
			for(int m=0; m<adjacent.xdim; m++){
				num+=dAij(adjacent, i, m)*dAij(adjacent, j, m);
				den+=max(dAij(adjacent, i, m), dAij(adjacent, j, m));
			}
			dAij(SimOut, i, j) = 2*num/den;
		}
	}
}

void HDIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	for(int i=0; i<SimOut.ydim; i++){
		for(int j=0; j<SimOut.xdim; j++){
			if(i==j){
				dAij(SimOut, i, j) = 1;
				continue;
			}
			if(j<i){
				dAij(SimOut, i, j) = dAij(SimOut, j, i);
				continue;
			}
			double intersection=0;
			double ki=0, kj=0;
			for(int m=0; m<adjacent.xdim; m++){
				intersection+=dAij(adjacent, i, m)*dAij(adjacent, j, m);
				ki+=dAij(adjacent, i, m);
				kj+=dAij(adjacent, j, m);
			}
			dAij(SimOut, i, j) = intersection/max(ki,kj);
		}
	}
}

void RAIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	for(int i=0; i<SimOut.ydim; i++){
		for(int j=0; j<SimOut.xdim; j++){
			if(i==j){
				dAij(SimOut, i, j) = 1;
				continue;
			}
			if(j<i){
				dAij(SimOut, i, j) = dAij(SimOut, j, i);
				continue;
			}
			dAij(SimOut, i, j) = 0;
			for(int m=0; m<adjacent.xdim; m++){
				if(dAij(adjacent, i, m)*dAij(adjacent, j, m)){
					int km=0;
					for(int k=0;k<adjacent.xdim; k++) 
						km+=dAij(adjacent,m,k);
					if(km) 
						dAij(SimOut,i,j)+=1.0/km;
				}
			}
		}
	}
}

void LPIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp, double eps)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	Matrix2D<double> A(adjacent.ydim, adjacent.xdim);
	adjacent.getSliceAsMatrix(0, A);
	Matrix2D<double> t=A*A;
	A=t+eps*A*t;

	SimOut=A;
}

void KatzIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp, double beta)
{
	SimOut.initZeros(SimIn);
	MultidimArray<double> adjacent(SimIn.ydim, SimIn.xdim);
	CreateAdjacentMatrix(SimIn, adjacent, k, cmp);

	Matrix2D<double> A(adjacent.ydim, adjacent.xdim), I(adjacent.ydim, adjacent.xdim);
	adjacent.getSliceAsMatrix(0, A);

	Matrix2D<double> U,V;
	Matrix1D<double> W;Matrix1D<int> index;
	A.eigs(U,W,V,index);
	beta=1/W(index(0))/2;

	A=beta*A;
	I.initIdentity();
	A=I-A;
	Matrix2D<double> Ainv;
	A.inv(Ainv);
	Ainv=Ainv-I;
	SimOut=Ainv;
}