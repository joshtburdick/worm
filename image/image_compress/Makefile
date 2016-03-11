
all : bwtCompress countRuns

bwtCompress : bwtCompress.cpp RunCompressorElias.h RunDecompressorElias.h EliasGammaCoder.h BWTEliasCompressor.h
	g++ -Wall -g -O6 -o bwtCompress bwtCompress.cpp /usr/local/lib/libsais.a

countRuns : countRuns.cpp CountRuns.h RunCounter.h
	g++ -Wall -g -o countRuns countRuns.cpp /usr/local/lib/libsais.a

