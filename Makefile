OUT    = ./bin/main
CC     = g++
INCDIR = ./include
OBJDIR = ./object
SRCDIR = .

SRCS   = $(shell find $(SRCDIR) -name "*.cpp")
OBJS   = $(SRCS:%.cpp=$(OBJDIR)/%.o)

FLAGS  = -fopenmp -Wall

obj:=main.o matrix.o

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
	@$(CC) $(FLAGS) -I$(INCDIR) -o $(OUT) $(OBJS)

$(OBJDIR)/%.o: %.cpp
	@echo "Compiling source: $<"
	@$(CC) $(FLAGS) -I$(INCDIR) -c $<  -o $@


.PHONY:clean

clean:
	rm -f $(OBJS) $(OUT)

