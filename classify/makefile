# program executable file name
PROG = classify

# list of source files
SRCS = $(wildcard *.cpp)
# list of generated object files
OBJS = $(SRCS:.cpp=.o)
# list of lib files
LIBS = -lgsl -lgslcblas -lm

# output directory
OUTDIR = bin
# boost directory
# BOOSTDIR = /usr/include/boost-1_33_1/
BOOSTDIR = ../boost_1_34_1/
# jsc directory
# JSCDIR = /home/jdu/programming/jdu_source_collection/
JSCDIR = ../jdu_source_collection/
# directory for INCLUDE searching 
# INCDIR = $(JSCDIR) $(BOOSTDIR)
INCDIR = $(JSCDIR) $(BOOSTDIR) ../cppunit/include/ ../gsl/include/
# directory for LIB searching
# LIBDIR =
LIBDIR = ../cppunit/lib/ ../gsl/lib/

# vpath settings
vpath %.h $(INCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OUTDIR)
vpath %   $(OUTDIR)

# target specific options
# -Wall			Enable all preprocessor warnings
# -w			Suppress warnings
# -g			Generate default debug format output
# -s --strip-all	Strip all symbols
# -mwindows		Create GUI application
# -mconsole		Create console application
# -Os			Optimise for space rather than speed
# -O[number]		Set optimization level to [number]
# -fconserve-space	Reduce size of object files
# -static		Do not link against shared libraries
# -mno-cygwin   build mingw app that do not depend on cygwin dll
# objdump -p a.exe | grep "DLL Name"
# TARGET_OPTION  = -Wall -O2
TARGET_OPTION  = -Wall -O2

# compile settings
CC	= g++
LD	= g++
CCFLAGS	= $(TARGET_OPTION) $(addprefix -I,$(INCDIR))
LDFLAGS	= $(TARGET_OPTION) $(addprefix -L,$(LIBDIR)) $(LIBS)

# command used to create directory
MKDIR	= mkdir
# command used to erase files
RM	= rm

# top-level rule, to compile everything
all: $(OUTDIR) $(PROG)

# rule to make the output directory if necessary
$(OUTDIR): 
	if [ ! -d "$(OUTDIR)" ]; then $(MKDIR) $(OUTDIR); fi

# rule to link the program
$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $(addprefix $(OUTDIR)/,$(OBJS)) -o $(OUTDIR)/$(PROG)

# meta rule for compiling any CPP source file
%.o: %.cpp
	$(CC) $(CCFLAGS) -c $< -o $(OUTDIR)/$@

# rule for cleaning re-compilable files
.PHONY: clean 
clean:
	$(RM) $(OUTDIR)/*
	$(RM) ./*.stackdump

# rule for generating the dependency file
.PHONY: depend 
depend:
	$(CC) $(CCFLAGS) -MM $(SRCS) > _dependency

sinclude _dependency
