original:
	g++ main.cpp -DORIGINAL -o original.out -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -I ../libs/opengm-master/src/external/TRWS-v1.3.src-patched/ -lexternal-library-qpbo-shared -lexternal-library-trws-shared  -Wfatal-errors
alp_exp:
	g++ main.cpp -DAE -o ae.out -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -I ../libs/opengm-master/src/external/TRWS-v1.3.src-patched/ -lexternal-library-qpbo-shared -lexternal-library-trws-shared  -Wfatal-errors
sym_ae:
	g++ main.cpp -DSYM_AE -o sym.out -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -I ../libs/opengm-master/src/external/TRWS-v1.3.src-patched/ -lexternal-library-qpbo-shared -lexternal-library-trws-shared  -Wfatal-errors
hybrid:
	g++ main.cpp -DHYBRID -o hyb.out -std=c++11 -lpng -I ../libs/opengm-master/src/external/QPBO-v1.3.src-patched/ -I ../libs/opengm-master/src/external/TRWS-v1.3.src-patched/ -lexternal-library-qpbo-shared -lexternal-library-trws-shared  -Wfatal-errors
all: alp_exp sym_ae
clean:
	rm *.out