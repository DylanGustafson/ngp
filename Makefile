.PHONY: ngp

CC = gcc
LINKGMP = -lgmp

ifeq ($(shell uname),Darwin)
	#Somewhat hacky way of getting the latest installed version of (actual) gcc
	CC = $(shell ls /opt/homebrew/bin | grep "^gcc-[0-9][0-9]" | tail -n 1)
	LINKGMP = -I /opt/homebrew/include -L /opt/homebrew/lib -lgmp
	MASM =
else
	MASM = -masm=intel
endif

ngp:
	rm -f ngp-bin psieve
	$(CC) -o psieve -O3 -march=native psieve32.c
	$(CC) -o ngp-bin -O3 -fopenmp -march=native $(MASM) ngp64.c $(LINKGMP)
