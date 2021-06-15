all: postProcess_v3.out

postProcess_v3.out: postProcess_v3.o collection1D.o data1D.o logger.o
	g++ -o postProcess_v3.out data1D.o collection1D.o logger.o postProcess_v3.o

postProcess_v3.o: postProcess_v3.cpp collection1D.h data1D.h logger.h
	g++ -c postProcess_v3.cpp 

collection1D.o: collection1D.cpp collection1D.h
	g++ -c collection1D.cpp

data1D.o: data1D.cpp data1D.h
	g++ -c data1D.cpp

logger.o: logger.cpp logger.h
	g++ -c logger.cpp

clean:
	rm -f *.o*
