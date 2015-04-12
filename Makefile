CC = gcc
CPP = g++

CFLAGS = -O3 -Iinclude
CPPFLAGS = $(CFLAGS)

LDFLAGS = -lgsl

OBJ = check_outliers.o gelman_rubin.o gen_CR.o \
			restore_state.o dream_initialize.o dream.o dream_pars.o

dream.a: $(OBJ)
	ar rcs dream.a $(OBJ)

clean:
	rm -rf dream.a
	rm -rf $(OBJ)

test_dream: test_dream.cpp dream.a
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

