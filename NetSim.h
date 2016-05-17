#ifndef __NETSIM_H__
#define __NETSIM_H__

#include <iostream>
#include <queue>
#include <vector>

#include <data/xmipp_image.h>

struct Elem{
	int index;
	double value;
	Elem(int i, double v){index=i;value=v;}
};

struct Cmp{
	bool reverse=false;
	Cmp(bool isreverse=false){reverse=isreverse;}
	bool operator()(Elem l, Elem r){
		if(reverse)
			return l.value>r.value;
		else
			return l.value<r.value;
	}	
};

void CreateAdjacentMatrix(MultidimArray<double> &SimIn, MultidimArray<double> &Out, int k, const Cmp cmp);
void CreateAdjacentMatrix(std::vector<double> &SimIn, std::vector<double> &Out, int k, const Cmp cmp);

//local Similarity measurement
void SNN(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp);
void JaccardIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp);
void SorensenIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp);
void HDIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp);
void RAIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp);
void LPIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOut, int k, const Cmp cmp, double eps);

//global similarity measurement
void KatzIndex(MultidimArray<double> &SimIn, MultidimArray<double> &SimOUt, int k, const Cmp cmp, double beta);

#endif