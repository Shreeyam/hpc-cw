EXE := ldc.out

SRC_DIR := src
OBJ_DIR := obj

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CC:=g++
CPPFLAGS := -I./include
CXXFLAGS := --std=c++11 -Wall -O2
LDFLAGS  := 
LDLIBS   := -llapack -lblas

DEBUGFLAGS:=
PRODFLAGS:=

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir $@

clean:
	-rm -f $(EXE)
	-rm -f $(OBJ_DIR)/*.o

debug: clean $(EXE)
	./$(EXE) $(DEBUGFLAGS)

production: clean $(EXE)
	./$(EXE) $(PRODFLAGS)