regular: *.c *.h
	gcc -fPIC -shared argsort.c -o libsort.so -lm
	mv -f libsort.so /usr/local/lib/
	gcc -fPIC -shared kdtree.c -o libkdtree.so -lm -lsort
	mv -f libkdtree.so /usr/local/lib/
	gcc -fPIC -shared energy.c -o libenergy.so -lm 
	mv -f libenergy.so /usr/local/lib/
	gcc -fPIC -shared kingmodel.c -o libking.so -lm 
	mv -f libking.so /usr/local/lib/
	gcc -fPIC -shared cluster_density.c -o libdensity.so -lm -lsort -lkdtree
	mv -f libdensity.so /usr/local/lib/
	cp -f *.h /usr/local/include/
local: *.c *.h
	gcc -fPIC -shared argsort.c -o libsort.so -lm
	gcc -fPIC -shared kdtree.c -o libkdtree.so -lm -lsort
	gcc -fPIC -shared energy.c -o libenergy.so -lm 
	gcc -fPIC -shared kingmodel.c -o libking.so -lm 
	gcc -fPIC -shared cluster_density.c -o libdensity.so -lm -lsort -lkdtree
