EXE := ldc.out

SRC_DIR := src
OBJ_DIR := obj

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CC:=g++
CPPFLAGS := -Iinclude
CFLAGS   := -std=c++1 -Wall -O2
LDFLAGS  := 
LDLIBS   := -llapack -lblas

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir $@

clean:
	-rm -f $(EXE)
	-rm -f $(OBJ_DIR)/*.o