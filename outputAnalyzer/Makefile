# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

all: analyzeOutput cleanup

analyzeOutput: analyzeOutput.o
	g++ -o analyzeOutput analyzeOutput.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

analyzeOutput.o: analyzeOutput.cpp
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o analyzeOutput.o analyzeOutput.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

cleanup:
	rm -rf analyzeOutput.o

clean:
	rm -rf analyzeOutput
	rm -rf analyzeOutput.o
