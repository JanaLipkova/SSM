CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O3 -std=c99 -Wall -D_XOPEN_SOURCE=700 -D_BSD_SOURCE

CFLAGS += `torc_cflags`

# GSL
#LDFLAGS += -lgsl -lgslcblas

# PRIMARY
LDFLAGS += `torc_libs` -lpthread -lm 


TARGETS += simulate_all

all: $(TARGETS)




simulate_all: propagation_tool.c  spawner.h  run_ssm.c  spawner.h
	$(CC) $(CFLAGS) -o simulate_all  propagation_tool.c $(LDFLAGS)





clean:
	rm -f simulate_all   *.o
