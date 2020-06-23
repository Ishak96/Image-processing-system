# name of output program
TARGET=object_tracking

CC=gcc

# projects dirs and configuration files
BINDIR=bin
SRCDIR=src
INCDIR=include
TESTDIR=test
OBJDIR=obj

CFLAGS=-Wall -I$(INCDIR) -w -g 
# getting source files and obj names
SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
TESTS    := $(wildcard $(TESTDIR)/*.c)
TESTBINS := $(TESTS:$(TESTDIR)/%.c=$(BINDIR)/%)
TESTOBJ  := $(TESTS:$(TESTDIR)/%.c=$(OBJDIR)/%.o)

# flags to be used
LFLAGS=-lm -std=c99 -lm

GTKFLAGS=`pkg-config --cflags gtk+-2.0 pkg-config --libs gtk+-2.0`

# make main executable
$(BINDIR)/$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CC) $^ $(LFLAGS) -o $@ $(GTKFLAGS)

# make obj files	$(CC) $^ $(LFLAGS) -o $@	$(CC) $^ $(LFLAGS) -o $@
$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@ $(GTKFLAGS)

# compile tests
test: $(TESTBINS)

$(TESTBINS): $(TESTOBJ) $(OBJECTS) | $(BINDIR)
	$(CC) $(filter-out $(OBJDIR)/main.o $(filter-out $(patsubst $(BINDIR)/%,$(OBJDIR)/%.o,$@), $(TESTOBJ)), $(OBJECTS) $(TESTOBJ)) -o $@ $(LFLAGS) $(GTKFLAGS)

$(TESTOBJ): $(TESTS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $(patsubst $(OBJDIR)/%.o,$(TESTDIR)/%.c,$@) -o $@ $(GTKFLAGS)

# Create the directories if they don't exist
$(OBJDIR):
	mkdir -p $@

$(BINDIR):
	mkdir -p $@


.PHONY: clean
clean:
	rm -f obj/* bin/*