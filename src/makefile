CC=g++-4.6
include ../makefile.inc
BIN=../bin

all: $(BIN)/ims-bicluster

ims-biclust: lib_json/json_reader.o lib_json/json_value.o lib_json/json_writer.o

$(BIN)/ims-bicluster: ims-bicluster.cpp 
	$(CC) $(CFLAGS) ims-bicluster.cpp $(LIB) -o $(BIN)/ims-bicluster

clean:
	rm -rf *o $(BIN)/ims-bicluster

