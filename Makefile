.PHONY: all clean

CC     := gcc
CFLAGS := -O3 -DHAVE_INLINE
LIBS   := -lgsl -lm -lblas
DEBUG  := # -ggdb -DDEBUG

SRCS := $(wildcard ./*.c)
HEADERS := $(wildcard ./*.h)

all: gmmem

clean:
	rm -f gmmem

gmmem: $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBS) $(DEBUG) -o $@ $(SRCS)
