include Make.inc

CPPFLAGS += -std=c++14 -g -m64 -O3
CPPFLAGS += -Iinclude -I.

LDFLAGS += -lgsl

SRC = $(wildcard *.cpp)

OBJ = $(SRC:%.cpp=build/%.o)

all: build/libdream.a

build/libdream.a: $(OBJ)
	@mkdir -p build
	ar rcs build/libdream.a $(OBJ)

clean:
	rm -rf build

build/test_dream: tests/test_dream.cpp build/libdream.a
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

##############################################################################

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

build/%.o: %.cpp
	@mkdir -p build
	$(CPP) $(CPPFLAGS) -o $@ -c $<

