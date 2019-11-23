DIST = README data doc src 
VERSION = 5-1

tgz :
	cd src; make clean
	tar zcvf mismatch$(VERSION).tgz README data src

zip :
	cd src; make clean
	zip -k -l -r mismatch$(VERSION).zip $(DIST) -x src/RCS/\* src/\*.o

clean :
	rm -rf *.zip *.tgz *~

