COMPILER = g++
BIN = bin
SRC = src
LIBS = -lm
CCFLAGS = -Wno-write-strings -O3 -mtune=native
OBJ = obj

all: MSTMain

MSTMain: Makefile ${SRC}/MSTMain.cpp ${SRC}/MSTMain.h $(OBJ)/MSTreeNode.o $(OBJ)/MSTree.o $(OBJ)/NormD.o
	${COMPILER} ${CCFLAGS} ${LIBS} -I /Users/WeiYang/software/gsl/gsl-1.16/include -L /Users/WeiYang/software/gsl/gsl-1.16/lib -lgsl -lgslcblas -o ${BIN}/MLnoPrior ${SRC}/MSTMain.cpp $(OBJ)/MSTree.o $(OBJ)/MSTreeNode.o $(OBJ)/NormD.o

$(OBJ)/NormD.o: Makefile ${SRC}/NormD.cpp ${SRC}/NormD.h
	${COMPILER} ${CCFLAGS} -c ${SRC}/NormD.cpp -o $(OBJ)/NormD.o

$(OBJ)/MSTreeNode.o: Makefile ${SRC}/MSTreeNode.cpp ${SRC}/MSTreeNode.h
	${COMPILER} ${CCFLAGS} -c ${SRC}/MSTreeNode.cpp -o $(OBJ)/MSTreeNode.o

$(OBJ)/MSTree.o: Makefile ${SRC}/MSTree.cpp ${SRC}/MSTree.h
	${COMPILER} ${CCFLAGS} -c ${SRC}/MSTree.cpp -o $(OBJ)/MSTree.o


clean:
	rm -rf obj/*.o
	rm ./bin/*

