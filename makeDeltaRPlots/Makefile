# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

all: makeDeltaRPlots cleanup

makeDeltaRPlots: makeDeltaRPlots.o
	g++ -o makeDeltaRPlots makeDeltaRPlots.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

makeDeltaRPlots.o: makeDeltaRPlots.cpp
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o makeDeltaRPlots.o makeDeltaRPlots.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

cleanup:
	rm -rf makeDeltaRPlots.o

clean:
	rm -rf makeDeltaRPlots
	rm -rf makeDeltaRPlots.o
