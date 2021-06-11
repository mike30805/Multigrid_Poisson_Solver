OUT    = ./bin/main
CC     = g++ -std=c++11
INCDIR = ./include
OBJDIR = ./object
SRCDIR = .

SRCS   = $(shell find $(SRCDIR) -name "*.cpp")
OBJS   = $(SRCS:%.cpp=$(OBJDIR)/%.o)

GSL_LIB = -L/home/return/gsl/lib -lgsl -lgslcblas
GSL_CFLAGS  = -I/home/return/gsl/include

FLAGS  = -fopenmp -Wall -O3

all: $(OUT)
	@echo "=========================================="
	@echo "Compiler          : $(CC)"
	@echo "Using flags       : $(FLAGS)"
	@echo "Include directory : $(INCDIR)"
	@echo "Object directory  : $(OBJDIR)"
	@echo "Excutable file    : $(OUT)"
	@echo "Successfully Done!"


$(OUT): $(OBJS)
	@echo "=========================================="
	@echo "Linking Executable $(OUT)"
	@$(CC) $(FLAGS) -I$(INCDIR) -o $(OUT) $(OBJS) $(GSL_LIB)

$(OBJDIR)/%.o: %.cpp
	@echo "Compiling source: $<"
	@$(CC) $(FLAGS) $(GSL_CFLAGS) -I$(INCDIR) -c $<  -o $@


.PHONY:clean

clean:
	rm -f $(OBJS) $(OUT)

