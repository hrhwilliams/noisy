SDIR = ./src
IDIR = ./incl
ODIR = ./obj

CPP = clang++
LIBS = -lz
CPPFLAGS := -std=c++17 -I$(IDIR)

_HPP = encoder.hpp noise.hpp
HPP = $(patsubst %,$(IDIR)/%,$(_HPP))

_OBJ = encoder.o noise.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp
	@mkdir -p $(@D)
	$(CPP) -c $< -o $@ $(CPPFLAGS)

noise: $(OBJ)
	$(CPP) -o $@ $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(ODIR)/*.o
