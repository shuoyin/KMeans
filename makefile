all: KMeansDivRA

KMeansDivRA: KMeans_main.o NetSim.o
	g++ -o KMeansDivRA -g NetSim.o KMeans_main.o -L/home/yinshuo/usr/xmipp/lib -L/usr/lib/openmpi/lib -lXmippData -lXmippClassif -lXmippDimred -lXmippExternal -lXmippParallel -lXmippRecons -ltiff -lmpi -lmpi_cxx

KMeans_main.o: KMeans.h KMeans.cpp NetSim.h
	g++ -g -I/home/yinshuo/usr/xmipp/libraries -I/home/yinshuo/usr/xmipp/externals -I/home/yinshuo/usr/xmipp -I/usr/include/openmpi -c KMeans.cpp -L/home/yinshuo/usr/xmipp/lib -lXmippData -lXmippClassif -ltiff -lXmippExternal -lXmippParallel -lXmippRecons -L/usr/lib/openmpi/lib -lmpi -lmpi_cxx -o KMeans_main.o

NetSim.o: NetSim.h NetSim.cpp
	g++ -g -I/home/yinshuo/usr/xmipp -I/home/yinshuo/usr/xmipp/libraries -I/home/yinshuo/us/xmipp/external -c NetSim.cpp -L/home/yinshuo/usr/xmipp/lib  -lXmippData -lXmippClassif -lXmippDimred -lXmippExternal -lXmippInterface -lXmippSqliteExt -lXmippParallel -lXmippRecons -lXmippJNI -ltiff -o NetSim.o

clean:
	rm *.o KMeansDivRA