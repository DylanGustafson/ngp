.PHONY: ngp

ngp:
	rm -f ngp ps
	gcc -O3 psieve32.c -o ps
	gcc -O3 -fopenmp ngp.c -lgmp -o ngp
