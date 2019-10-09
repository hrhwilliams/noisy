SDIR = ./src
IDIR = ./incl
ODIR = ./obj

CPP = g++
LIBS = -lz
CPPFLAGS := -std=c++17 -I$(IDIR) -g

_HPP = encoder.hpp noise.hpp mask.hpp
HPP = $(patsubst %,$(IDIR)/%,$(_HPP))

_OBJ = encoder.o noise.o mask.o noise_functions.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp
	@mkdir -p $(@D)
	$(CPP) -c $< -o $@ $(CPPFLAGS)

noise: $(OBJ)
	$(CPP) -o $@ $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(ODIR)/*.o
