#include "KMeans.h"

using namespace std;

double fitTwoImages(MultidimArray<double> &fixed,MultidimArray<double> &I, double sigma, bool reverse)
{
	if (reverse)
    {
        I.selfReverseX();
        I.setXmippOrigin();
    }

    Matrix2D<double> M;

    alignImages(fixed, I, M);

    // Compute the correntropy
    double corr=0.0;
    corr = correntropy(fixed, I, sigma);

    return corr;
}

double fit(MultidimArray<double> &fixed, MultidimArray<double> &I, double sigma)
{
	MultidimArray<double> Idirect = I;
	double corrd = fitTwoImages(fixed, Idirect, sigma);

	MultidimArray<double> Imirror=I;
	double corri = fitTwoImages(fixed, Imirror, sigma, true);

	if(corrd > corri){
		I=Idirect;
		return corrd;
	}
	else{
		I=Imirror;
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
    I().statisticsAdjust(0, 1);
}

void printCenterInfo(KMeans &km, string name)
{
	for(int i=0; i<km.K; i++){
		Cluster &c = km.centers[i];
		Image<double> out=c.thisP;
		out.write(name+"_centers.stk", i, true, WRITE_APPEND);
		cout<<"size of Cluster "<< i<<": "<<c.belongings.size()<<endl;
		cout<<"it contains: ";
		for(int i=0; i<c.belongings.size(); i++)
			cout<<c.belongings[i]<<" ";
		cout<<endl<<endl;
	}
}

//member functions of class Program

Program::Program(string s)
	:fin(s)
{
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

    corrM.resize(Nimgs,Nimgs);
    for(int i=0; i<Nimgs; i++){
        Image<double> img1;
        readImage(md, img1, i+1, false);
        for(int j=0; j<Nimgs; j++){
            if(i==j){
                dAij(corrM, i, j)=1;
                continue;
            }
            if(j<i){
                dAij(corrM, i, j) = dAij(corrM, j, i);
                continue;
            }
            Image<double> img2;
            readImage(md, img2, j+1, false);
            dAij(corrM, i, j) = correntropy(img1.data, img2.data, sigma);
        }
    }
    MultidimArray<double> tmp=corrM;
    RAIndex(tmp, corrM, 80, Cmp(false));
    CreateAdjacentMatrix(corrM, adjancent, 80, Cmp(false));

}

bool Program::run()
{
	KMeans km(3,30, fin, Nimgs, Xdim, Ydim, sigma);
	km.initialize(adjancent);
	bool r=true;
	km.runDivisive(4,adjancent);
	cout<<"run done"<<endl;
	return r;
}

//member functions for KMeans
void KMeans::createCluster(Cluster &c)
{
	c.thisP.initZeros(Ydim, Xdim);
	c.nextP.initZeros(Ydim, Xdim);
	c.thisP.setXmippOrigin();
	c.nextP.setXmippOrigin();
	c.corrV.resize(Nimgs, -1);
	c.belongings.clear();
	c.newbelongings.clear();

}

void KMeans::initialize(MultidimArray<double> &adjancent)
{
	for(int i=0; i<K; i++){
		Cluster c;
		centers.push_back(c);
		createCluster(centers[i]);
	}

	for(int i=0; i<Nimgs; i++){
		Image<double> I;
		readImage(md, I, i+1, false);
		centers[i%K].nextP+=I();
		centers[i%K].newbelongings.push_back(i);
		newclusters.push_back(i%K);
	}

	for(int i=0; i<K; i++)
		centers[i].update(md, sigma, adjancent);
}

void KMeans::splitCluster(Cluster &node, Cluster &node1, Cluster &node2, MultidimArray<double> &adjancent)
{
	//clear two new clusters
	createCluster(node1);createCluster(node2);

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
	node1.update(md, sigma, adjancent);
	node2.update(md, sigma, adjancent);
}

bool KMeans::run(MultidimArray<double> &adjancent)
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
				double corr=centers[c].computeNetSim(i, adjancent);
				if(corr>bestcorr){
					bestcorr=corr;
					bestcluster=c;
				}
			}
			fit(centers[bestcluster].thisP, I(), sigma);
			centers[bestcluster].nextP+=I();
			centers[bestcluster].newbelongings.push_back(i);
			newclusters[i]=bestcluster;
		}
		for(int i=0; i<K; i++) centers[i].update(md, sigma, adjancent);
		int changes=0;
		for(int i=0; i<Nimgs; i++)
			changes+=(oldclusters[i]==newclusters[i] ? 0 : 1);
		cout<<"in iter "<<ite<<", changes "<<changes<<endl;
		printCenterInfo(*this, name.removeLastExtension().getString());
		cout<<endl;

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
			splitCluster(centers[largest], centers[K], centers[K+1], adjancent);
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

void KMeans::runDivisive(int finalN,MultidimArray<double> &adjancent)
{
	while(true){
		run(adjancent);
		printCenterInfo(*this, name.removeLastExtension().getString());
		if(K>=finalN)
			break;
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
		splitCluster(centers[largest], centers[K], centers[K+1], adjancent);
		centers.erase(centers.begin()+largest);
		K++;
		for(int i=0; i<K; i++){
			for(int j=0; j<centers[i].belongings.size(); j++)
				newclusters[centers[i].belongings[j]]=i;
		}
	}
}

//member function for Cluster
void Cluster::update(MetaData &md, double sigma, MultidimArray<double> &adjancent)
{
	belongings=newbelongings;
	newbelongings.clear();

	thisP = nextP/belongings.size();
	nextP.initZeros();

	int n = md.size();
	for(int i=0; i<n; i++){
		Image<double> I;
		readImage(md, I, i+1, false);
		corrV[i]=fit(thisP, I(), sigma);
	}
	vector<double> tmp=corrV;
	for(int i=0; i<n; i++)
		tmp[i]=computeNetSim(i, adjancent);
	corrV=tmp;
}

double Cluster::computeNetSim(int index, MultidimArray<double> &adjancent)
{
	int sn=80;
	std::vector<double> v;
	CreateAdjacentMatrix(corrV, v, sn, Cmp(false));
	double sim=0;
	for(int i=0; i<v.size(); i++){
		// sim+=v[i]*dAij(adjancent, index, i);
		if(dAij(adjancent, index, i)*v[i]){
			int km=0;
			for(int k=0;k<adjancent.xdim; k++) 
				km+=dAij(adjancent,i,k);
			if(km) 
				sim+=1.0/km;
		}
	}
	// sim/=(2*sn-sim);
	return sim;
}


int main(int argc, char**argv)
{
	cout<<"Input image: "<<argv[1]<<endl;
	Program p(argv[1]);
	p.run();
	return 0;
}