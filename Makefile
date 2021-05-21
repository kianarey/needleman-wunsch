NVCC        = nvcc

NVCC_FLAGS  = -I/usr/local/cuda/include -gencode=arch=compute_60,code=\"sm_60\"
ifdef dbg
	NVCC_FLAGS  += -g -G
else
	NVCC_FLAGS  += -O2
endif

LD_FLAGS    = -lcudart -L/usr/local/cuda/lib64
EXE	        = needleman_wunsch
OBJ	        = Needleman_Wunsch_CPU.o main.o util.o

default: $(EXE)

main.o: main.cu util.h util.cpp NW_Single_Kernel.cu
	$(NVCC) -c -o $@ main.cu $(NVCC_FLAGS) 	

Needleman_Wunsch.o: Needleman_Wunsch_CPU.cpp Needleman_Wunsch_CPU.h
	$(NVCC) -c -o $@ Needleman_Wunsch_CPU.cpp $(NVCC_FLAGS)

util.o: util.cpp util.h
	$(NVCC) -c -o $@ util.cpp $(NVCC_FLAGS) 

$(EXE): $(OBJ)
	$(NVCC) $(OBJ) -o $(EXE) $(LD_FLAGS) $(NVCC_FLAGS)

clean:
	rm -rf *.o $(EXE)
