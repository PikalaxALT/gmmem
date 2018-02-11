.PHONY: all clean

CC     := gcc
CFLAGS := -O3
LIBS   := -lgsl -lm -lblas
DEBUG  := -ggdb -DDEBUG -DHAVE_INLINE

SRCS := $(wildcard ./*.c)

all: gmmem

clean:
	rm -f gmmem

gmmem: $(SRCS)
	$(CC) $(CFLAGS) $(LIBS) $(DEBUG) -o $@ $^
