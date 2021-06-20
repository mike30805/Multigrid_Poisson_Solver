OUT    = ./bin/main
CC     = g++ -std=c++11 #CPU
NVCC     = nvcc -std=c++11 #GPU
INCDIR = ./include
OBJDIR = ./object
SRCDIR = .
SIMU_OPTION =
# simulation options
#######################################################################################################

# (a) GPU acceleration
# --> must turn off OMP
#SIMU_OPTION += -DGPU

# (b) OMP acceleration
# --> must turn off GPU
SIMU_OPTION += -DOMP

SRCS  := $(shell find $(SRCDIR) -name "*.cpp")
OBJS   = $(SRCS:%.cpp=$(OBJDIR)/%.o)
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
SRCS  := $(filter-out ./matrix.cpp, $(SRCS))
OBJS  := $(filter-out ./matrix.o, $(OBJS))
endif


ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
SRCS_GPU   = $(shell find $(SRCDIR) -name "*.cu")
OBJS_GPU   = $(SRCS_GPU:%.cu=$(OBJDIR)/%.o)
OBJ_GPU_LINK := $(OBJDIR)/gpu_link.o
endif

GSL_LIB = -L/software/gsl/default/lib -lgsl -lgslcblas
GSL_CFLAGS  = -I/software/gsl/default/include



ifeq "$(filter -DOMP, $(SIMU_OPTION))" "-DOMP"
FLAGS  = -fopenmp -Wall #CPU
else
FLAGS  = -Wall #CPU
endif

ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
NVCC_FLAGS  =   #GPU
endif

all: $(OUT)
	@echo "=========================================="
	@echo "Compiler          : $(CC)"
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
	@echo "GPU Compiler       : $(NVCC)"
endif
	@echo "Using flags       : $(FLAGS)"
	@echo "Include directory : $(INCDIR)"
	@echo "Object directory  : $(OBJDIR)"
	@echo "Excutable file    : $(OUT)"
	@echo "Successfully Done!"



$(OUT): $(OBJS) $(OBJS_GPU)
	@echo "=========================================="
	@echo "Linking Executable $(OUT)"
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
	@echo "Linking GPU codes"
	@$(NVCC) -o $(OBJ_GPU_LINK) $(OBJS_GPU) -dlink
endif
	@echo "Linking CPU codes"
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
	@$(NVCC) -o $@ $^ $(OBJ_GPU_LINK) $(FLAGS_GPU) $(GSL_LIB)
else
	@$(CC) $(FLAGS) -I$(INCDIR) -o $(OUT) $(OBJS)  $(GSL_LIB)
endif



$(OBJDIR)/%.o: %.cpp 
	@echo "Compiling source: $<"
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
	@if [ $< !=  "matrix.cpp" ]; \
    then $(CC) $(FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@ ;\
	else $(NVCC) $(NVCC_FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c matrix.cu  -o $@;\
	fi
else
	@$(CC) $(FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@
endif

$(OBJDIR)/%.o: %.cu
	@echo "Compiling source: $<"
	@$(NVCC) $(NVCC_FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@



.PHONY:clean

clean:
	rm -f $(OBJDIR)/*.o

