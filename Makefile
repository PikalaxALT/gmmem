.PHONY: all clean

CC     := gcc
CFLAGS := -O3 -DHAVE_INLINE
LIBS   := -lgsl -lm -lblas
DEBUG  := # -ggdb -DDEBUG

SRCS := $(wildcard ./*.c)

all: gmmem

clean:
	rm -f gmmem

gmmem: $(SRCS)
	$(CC) $(CFLAGS) $(LIBS) $(DEBUG) -o $@ $^
