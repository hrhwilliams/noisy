CPP := clang++
SDIR = ./src
IDIR = ./incl
CPPFLAGS := -std=c++17 -lz -I$(IDIR)

all:
	$(CPP) $(SDIR)/main.cpp $(SDIR)/noise.cpp $(SDIR)/encoder.cpp $(CPPFLAGS) -o Noise
