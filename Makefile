.PHONY : all clean

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -funroll-loops  -DNDEBUG -DBOOST
DEBUGFLAGS = -std=c++17 -g -pg -Wall -funroll-loops -DNDEBUG 
INCLUDE_TARGETS=-Iinclude/ -I./ -I/opt/homebrew/opt/boost/include/
LIBS=-L/software/boost/1.78.0/b1/lib/ -lboost_program_options
# GSET=/scratch/mhuang_lab/ISING_MACHINES/GSET/tiny_new/G281
# BEST=28
# GSET=/scratch/mhuang_lab/ISING_MACHINES/GSET/K_graphs/K16384
# BEST=790000
TARGET=pt_icm
GRAPH=$(GSET)/set/G001
BEST=11624

# GSET=/scratch/mhuang_lab/ISING_MACHINES/GSET/limit_int_sparse/G001
# BEST=125

# GSET=/scratch/mhuang_lab/ISING_MACHINES/GSET/SATLIB_GRAPHS/SAT_TO_ISING/uf20-91/uf20-01
# BEST=47

INIT=/scratch/mhuang_lab/ISING_MACHINES/INITS/RND/000.txt

# GSET=/scratch/mhuang_lab/ashar36/chimera_1_2_4
# BEST=0
# INIT=/scratch/mhuang_lab/ashar36/init_chim.txt

SRC := $(wildcard src/*.cc)
OBJ := $(patsubst src/%.cc,bin/%.o, $(SRC))

all : $(TARGET)

bin/*.o  : include/*.hh

$(TARGET) : $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE_TARGETS) -o $@ $^ $(LIBS)

bin/%.o : src/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE_TARGETS) $(MACROS) -o $@ -c $< $(LIBS)
run: $(TARGET)
	./$(TARGET) $(GRAPH)
debug: CXXFLAGS=$(DEBUGFLAGS)
debug: $(OBJ)
	$(CXX) $(INCLUDE_TARGETS) $(CXXFLAGS) -o $@ $^  $(LIBS)

clean :
	rm -rf $(OBJ) $(TARGET) debug
