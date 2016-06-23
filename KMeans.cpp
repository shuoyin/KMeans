#include "KMeans.h"
#include <sstream>

using namespace std;
using namespace MPI;

string pro;

double fitTwoImages(MultidimArray<double> &fixed,MultidimArray<double> &I, double sigma, Matrix2D<double> &M, bool reverse)
{
	if (reverse)
    {
        I.selfReverseX();
        I.setXmippOrigin();
    }

    alignImages(fixed, I, M);

    // Compute the correntropy
    double corr=0.0;
    corr = correntropy(fixed, I, sqrt(2)*sigma);

    return corr;
}

double fit(MultidimArray<double> &fixed, MultidimArray<double> &I, double sigma, Matrix2D<double> &M)
{
	MultidimArray<double> Idirect = I;
	Matrix2D<double> md;
	double corrd = fitTwoImages(fixed, Idirect, sigma, md);

	MultidimArray<double> Imirror=I;
	Matrix2D<double> mi;
	// double corri=-100;
	double corri = fitTwoImages(fixed, Imirror, sigma, mi, true);

	if(corrd > corri){
		I=Idirect;
		M=md;
		return corrd;
	}
	else{
		I=Imirror;
		M=mi;
		return corri;
	}
}

void readImage(MetaData &md, Image<double> &I, size_t objId, bool applyGeo)
{
    if (applyGeo)
        I.readApplyGeo(md, objId);
    else
    {
        FileName fnImg;
        md.getValue(MDL_IMAGE, fnImg, objId);
        I.read(fnImg);
    }
    I().setXmippOrigin();
    // realGaussianFilter(I(), 1);
    I().statisticsAdjust(0, 1);
}

void printCenterInfo(KMeans &km, string name)
{
	for(int i=0; i<km.K; i++){
		Cluster &c = km.centers[i];
		Image<double> out=c.thisP;
		out.write(name+"_"+pro+"_centers.stk", i, true, WRITE_APPEND);
		cout<<"size of Cluster "<< i<<": "<<c.belongings.size()<<endl;
		cout<<"it contains: ";
		for(int i=0; i<c.belongings.size(); i++)
			cout<<c.belongings[i]<<" ";
		cout<<endl<<endl;
	}
}

//member functions of class Program

Program::Program(string s, int k0, int kn)
	:fin(s), K0(k0), KN(kn)
{
	time_t start,end;
    start=time(NULL);

	md.read(fin);

	Nimgs = md.size();
	size_t Zdim, Ndim;
    getImageSize(md, Xdim, Ydim, Zdim, Ndim);

	cout<<"Computing sigma..."<<endl;
	MultidimArray<int> mask(Ydim,Xdim);
	mask.setXmippOrigin();
    BinaryCircularMask(mask, Xdim/2, INNER_MASK);

    for(int i=1; i<=Nimgs; i++){
    	Image<double> I;
    	readImage(md, I, i, true);
    	double avg, stddev;
        I().computeAvgStdev_within_binary_mask(mask, avg, stddev);
        sigma += stddev;
    }
    sigma *= sqrt(2.0) / Nimgs;

    cout<<"Computing corrM..."<<endl;

    // corrM.resize(Nimgs,Nimgs);
    // for(int i=0; i<Nimgs; i++){
    // 	if(i%numprocs!=id)
    // 		continue;
    //     Image<double> img1;
    //     readImage(md, img1, i+1, false);
    //     for(int j=0; j<Nimgs; j++){
    //         if(i==j){
    //             dAij(corrM, i, j)=1;
    //             continue;
    //         }
    //         // if(j<i){
    //         //     dAij(corrM, i, j) = dAij(corrM, j, i);
    //         //     continue;
    //         // }
    //         Image<double> img2;
    //         readImage(md, img2, j+1, false);
    //         dAij(corrM, i, j) = correntropy(img1.data, img2.data, sqrt(2)*sigma);
    //     }
    // }
    // comm.Barrier();
    // for(int i=0; i<Nimgs; i++){
    // 	int source = i%numprocs;
    // 	comm.Bcast(&dAij(corrM, i, 0), corrM.xdim, MPI_DOUBLE, source);
    	
    // }
    // // MultidimArray<double> tmp=corrM;
    // // RAIndex(tmp, corrM, 80, Cmp(false));
    // Image<double> out=corrM;out.write("corr.stk");
    // // CreateAdjacentMatrix(corrM, adjancent, 80, Cmp(false));
    Image<double> I;
    I.read("corr.stk");
    corrM=I();

    end=time(NULL);
    cout<<"cost "<<difftime(end, start)<<" seconds"<<endl;
}

bool Program::run()
{
	KMeans km(K0,20, fin, Nimgs, Xdim, Ydim, sigma);
	km.initialize(corrM);
	bool r=true;
	km.runDivisive(KN,corrM);
	cout<<"run done"<<endl;
	return r;
}

//member functions for KMeans
void KMeans::createCluster(Cluster &c, MultidimArray<double> &corrM)
{
	c.thisP.initZeros(Ydim, Xdim);
	c.nextP.initZeros(Ydim, Xdim);
	c.thisP.setXmippOrigin();
	c.nextP.setXmippOrigin();
	c.corrV.resize(Nimgs, -1);
	c.belongings.clear();
	c.newbelongings.clear();

	int dim=Nimgs+1;
	c.simM.initZeros(dim,dim);c.netSimM.initZeros(dim,dim);
	for(int i=0; i<dim; i++){
		if(i<dim-1){
			for(int j=0; j<dim; j++){
				if(j==dim-1)
					dAij(c.simM, i, j)=-1;
				else
					dAij(c.simM, i, j)=dAij(corrM, i, j);
			}
		}
		else{
			for(int j=0; j<dim; j++) dAij(c.simM, i, j)=-1;
		}
	}
}

void KMeans::initialize(MultidimArray<double> &corrM)
{
	for(int i=0; i<K; i++){
		Cluster c;
		centers.push_back(c);
		createCluster(centers[i], corrM);
	}

	for(int i=0; i<Nimgs; i++){
		Image<double> I;
		readImage(md, I, i+1, false);
		centers[i%K].nextP+=I();
		centers[i%K].newbelongings.push_back(i);
		newclusters.push_back(i%K);
	}

	for(int i=0; i<K; i++)
		centers[i].update(md, sigma);
}

void KMeans::splitCluster(Cluster &node, Cluster &node1, Cluster &node2, MultidimArray<double> &corrM)
{
	//clear two new clusters
	createCluster(node1, corrM);createCluster(node2, corrM);

	int n=node.belongings.size();

	priority_queue<Elem, vector<Elem>, Cmp> q(Cmp(false));
	for(int i=0; i<n; i++){
		int index=node.belongings[i];
		q.push(Elem(index, node.corrV[index]));
	}

	for(int i=0; i<n/2; i++){
		int index = q.top().index;
		q.pop();
		Image<double> I;
		readImage(md, I, index+1, false);
		node1.nextP+=I();
		node1.newbelongings.push_back(index);
	}

	for(int i=n/2; i<n; i++){
		int index = q.top().index;
		q.pop();
		Image<double> I;
		readImage(md, I, index+1, false);
		node2.nextP+=I();
		node2.newbelongings.push_back(index);
	}
	node1.update(md, sigma);
	node2.update(md, sigma);
}

bool KMeans::run(MultidimArray<double> &corrM)
{
	for(int ite=0; ite<iter; ite++){
		oldclusters=newclusters;
		for(int i=0; i<newclusters.size(); i++) newclusters[i]=-1;
		for(int i=0; i<Nimgs; i++){
			Image<double> I;
			readImage(md, I, i+1, false);
			double bestcorr=-1;
			int bestcluster=-1;
			for(int c=0; c<K; c++){
				double corr=centers[c].getNetSim(i);
				if(corr>bestcorr){
					bestcorr=corr;
					bestcluster=c;
				}
			}
			Matrix2D<double> M;
			fit(centers[bestcluster].thisP, I(), sigma, M);
			centers[bestcluster].nextP+=I();
			centers[bestcluster].newbelongings.push_back(i);
			newclusters[i]=bestcluster;
		}

		for(int i=0; i<K; i++) centers[i].update(md, sigma);
		int changes=0;
		for(int i=0; i<Nimgs; i++)
			changes+=(oldclusters[i]==newclusters[i] ? 0 : 1);
		if(id==0){
			cout<<"in iter "<<ite<<", changes "<<changes<<endl;
			printCenterInfo(*this, name.removeLastExtension().getString());
			cout<<endl;
		}

		int largest=-1, largestsize=0, smallest=-1, smallsize=Nimgs+1;
		for(int i=0; i<K; i++)
		{
			if(centers[i].belongings.size()>largestsize){
				largestsize=centers[i].belongings.size();
				largest=i;
			}
			if(centers[i].belongings.size()<smallsize){
				smallsize=centers[i].belongings.size();
				smallest=i;
			}
		}
		if(smallsize<20){
			Cluster c;
			centers.push_back(c);centers.push_back(c);
			splitCluster(centers[largest], centers[K], centers[K+1], corrM);
			if(largest>smallest){
				centers.erase(centers.begin()+largest);
				centers.erase(centers.begin()+smallest);
			}
			else{
				centers.erase(centers.begin()+smallest);
				centers.erase(centers.begin()+largest);
			}
		}
		if(changes<5 && smallsize>20)
			break;
	}
	return true;
}

void KMeans::runDivisive(int finalN,MultidimArray<double> &corrM)
{
	while(true){
		time_t start, end;
		start=time(NULL);
		run(corrM);
		end=time(NULL);
		cout<<"One KMeans cost "<<difftime(end, start)<<" seconds"<<endl;
		if (id==0) printCenterInfo(*this, name.removeLastExtension().getString());
		if(K>=finalN)
			break;
		int tosplit = min(K, finalN-K);
		for(int j=0; j<tosplit; j++){
			int largest=-1, largestsize=0;
			for(int i=0; i<K; i++)
			{
				if(centers[i].belongings.size()>largestsize){
					largestsize=centers[i].belongings.size();
					largest=i;
				}
			}
			Cluster c;
			centers.push_back(c);centers.push_back(c);
			splitCluster(centers[largest], centers[K], centers[K+1], corrM);
			centers.erase(centers.begin()+largest);
			K++;
			for(int i=0; i<K; i++){
				for(int j=0; j<centers[i].belongings.size(); j++)
					newclusters[centers[i].belongings[j]]=i;
			}
		}
	}
	writeResult();
}

void KMeans::writeResult()
{
	if(id!=0) return ;
	printCenterInfo(*this, name.removeLastExtension().getString()+"final");
	MetaData out;
	FileName fn=name.removeLastExtension().getString()+"_"+pro+"results.xmd";
	FileName classname=name.removeLastExtension().getString()+"_"+pro+"final_centers.stk";
	for(int i=0; i<K; i++){
		MDRow row;
		FileName fnclass;
		fnclass.compose(i+1, classname);
		row.setValue(MDL_REF, i+1);
		row.setValue(MDL_IMAGE, fnclass);
		row.setValue(MDL_CLASS_COUNT, centers[i].belongings.size());
		out.addRow(row);
	}
	out.write("classes@"+fn);

	for(int i=0; i<K; i++){
		MetaData oneclass;
		vector<int> &bel=centers[i].belongings;
		for(int j=0; j<bel.size(); j++){
			MDRow tmp;
			md.getRow(tmp, bel[j]+1);
			Matrix2D<double> M;
			Image<double> I;
			readImage(md, I, bel[j]+1, false);
			fit(centers[i].thisP, I(), sigma, M);
			bool flip;
			double shiftx, shifty, psi, scale;
			transformationMatrix2Parameters2D(M, flip, scale, shiftx, shifty, psi);
			tmp.setValue(MDL_FLIP, flip);
            tmp.setValue(MDL_SHIFT_X, shiftx);
            tmp.setValue(MDL_SHIFT_Y, shifty);
            tmp.setValue(MDL_ANGLE_PSI, psi);
            oneclass.addRow(tmp);
		}
		oneclass.write(formatString("class%06d_images@%s",i+1,fn.c_str()),MD_APPEND);
	}
}

//member function for Cluster
void Cluster::update(MetaData &md, double sigma)
{
	belongings=newbelongings;
	newbelongings.clear();

	thisP = nextP/belongings.size();
	nextP.initZeros();

	int n = md.size();
	for(int i=0; i<n; i++){
		if(i%numprocs!=id)
			continue;
		Image<double> I;
		readImage(md, I, i+1, false);
		Matrix2D<double> M;
		corrV[i]=fit(thisP, I(), sigma, M);
	}
	comm.Barrier();
	for(int i=0; i<n; i++){
		int source = i%numprocs;
		comm.Bcast(&corrV[i], 1, MPI_DOUBLE, source);
	}
	for(int i=0; i<n; i++) {
		dAij(simM, n, i)=corrV[i];
		dAij(simM, i, n)=corrV[i];
	}
	dAij(simM, n, n)=1;
	
	int sn=250;
	SorensenIndex(simM, netSimM, sn, Cmp(false));
	// MultidimArray<double> tmp=netSimM;
	// KatzIndex(tmp, netSimM, sn, Cmp(false), 0.01);
}

double Cluster::getNetSim(int index)
{
	// int sn=80;
	// std::vector<double> v;
	// CreateAdjacentMatrix(corrV, v, sn, Cmp(false));
	// double sim=0;
	// for(int i=0; i<v.size(); i++){
	// 	sim+=v[i]*dAij(adjancent, index, i);
	// }
	// sim/=(2*sn-sim);
	// return sim;
	return dAij(netSimM, netSimM.xdim-1, index);
}


int main(int argc, char**argv)
{
	time_t start, end;
	start=time(NULL);
	Init(argc, argv);
	id=comm.Get_rank();
	numprocs=comm.Get_size();
	cout<<"Input image: "<<argv[1]<<endl;
	stringstream c1(argv[2]), c2(argv[3]);
	int k0,kn;
	c1>>k0;c2>>kn;
	pro=argv[0];
	Program p(argv[1], k0, kn);
	p.run();
	Finalize();
	end=time(NULL);
	cout<<"all process cost "<<difftime(end, start)<<" seconds"<<endl;
	return 0;
}