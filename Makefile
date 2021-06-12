OUT    = ./bin/main
CC     = g++ -std=c++11 #CPU
NVCC     = nvcc -std=c++11 #GPU
INCDIR = ./include
OBJDIR = ./object
SRCDIR = .

SRCS   = $(shell find $(SRCDIR) -name "*.cpp")
OBJS   = $(SRCS:%.cpp=$(OBJDIR)/%.o)

SRCS_GPU   = $(shell find $(SRCDIR) -name "*.cu")
OBJS_GPU   = $(SRCS_GPU:%.cu=$(OBJDIR)/%.o)
OBJ_GPU_LINK := $(OBJDIR)/gpu_link.o

GSL_LIB = -L/software/gsl/default/lib -lgsl -lgslcblas
GSL_CFLAGS  = -I/software/gsl/default/include

FLAGS  = -fopenmp -Wall -O3 #CPU
NVCC_FLAGS  =  -O3 #GPU

all: $(OUT)
	@echo "=========================================="
	@echo "Compiler          : $(CC)"
	@echo "Using flags       : $(FLAGS)"
	@echo "Include directory : $(INCDIR)"
	@echo "Object directory  : $(OBJDIR)"
	@echo "Excutable file    : $(OUT)"
	@echo "Successfully Done!"


$(OUT): $(OBJS) $(OBJS_GPU)
	@echo "=========================================="
	@echo "Linking Executable $(OUT)"
	@$(NVCC) -o $(OBJ_GPU_LINK) $(OBJS_GPU) -dlink

	@$(CC) $(FLAGS) -I$(INCDIR) -o $(OUT) $(OBJS) $(OBJ_GPU_LINK)  $(GSL_LIB)

$(OBJDIR)/%.o: %.cpp
	@echo "Compiling source: $<"
	@$(CC) $(FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@

$(OBJDIR)/%.o: %.cu
	@echo "Compiling source: $<"
	@$(NVCC) $(NVCC_FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@


.PHONY:clean

clean:
	rm -f $(OBJS) $(OUT)

