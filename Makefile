.PHONY: ngp

ngp:
	rm -f ngp ps
	gcc -O3 -march=native psieve32.c -o ps
	gcc -O3 -fopenmp -march=native ngp64.c -lgmp -o ngp
