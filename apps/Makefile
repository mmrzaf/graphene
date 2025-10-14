CC = gcc
CFLAGS = -Wall -Wextra -O2 -I../include
LDFLAGS = -L../lib -lgraphene
TARGETS = example_app gwo_im optim_compare

all: $(TARGETS)

example_app: example_app.c
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

gwo_im: gwo_im.c
	$(CC) $(CFLAGS) $< $(LDFLAGS) -lm -o $@

optim_compare: optim_compare.c
	$(CC) $(CFLAGS) $< $(LDFLAGS) -lm -o $@

clean:
	rm -f $(TARGETS) gwo_im optim_compare

