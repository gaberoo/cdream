CC = gcc
CPP = g++

CFLAGS = -m64 -O3 -Iinclude
CPPFLAGS = $(CFLAGS)

LDFLAGS = -lgsl

OBJ = check_outliers.o gelman_rubin.o gen_CR.o \
			restore_state.o dream_initialize.o dream.o \
			dream_pars.o

libdream.a: $(OBJ)
	ar rcs libdream.a $(OBJ)

clean:
	rm -rf libdream.a
	rm -rf $(OBJ)

test_dream: test_dream.cpp libdream.a
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

