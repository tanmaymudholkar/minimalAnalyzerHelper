# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

all: makeTriggerEfficiencyPlots cleanup

makeTriggerEfficiencyPlots: makeTriggerEfficiencyPlots.o
	g++ -o makeTriggerEfficiencyPlots makeTriggerEfficiencyPlots.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

makeTriggerEfficiencyPlots.o: makeTriggerEfficiencyPlots.cpp
	g++ -g -c -std=c++1z -Wall -Wextra -Werror -pedantic -pedantic-errors -fPIC -O3 -o makeTriggerEfficiencyPlots.o makeTriggerEfficiencyPlots.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

cleanup:
	rm -rf makeTriggerEfficiencyPlots.o

clean:
	rm -rf makeTriggerEfficiencyPlots
	rm -rf makeTriggerEfficiencyPlots.o
