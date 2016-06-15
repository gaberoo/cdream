include Make.inc

CPPFLAGS += -std=c++14 -g -m64 -O3
CPPFLAGS += -Iinclude

LDFLAGS = -lgsl

SRC = $(wildcard *.cpp)

OBJ = $(SRC:%.cpp=build/%.o)

all: build/libdream.a

build/libdream.a: $(OBJ)
	@mkdir -p build
	ar rcs libdream.a $(OBJ)

clean:
	rm -rf build

test_dream: test_dream.cpp libdream.a
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

##############################################################################

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

build/%.o: %.cpp
	@mkdir -p build
	$(CPP) $(CPPFLAGS) -o $@ -c $<

