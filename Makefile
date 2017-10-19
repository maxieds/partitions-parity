#### Makefile 
#### Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
#### Updated: 2017.10.18 

CC = g++
EXE = gencerts
CFLAGS = -g -lm -std=c++11 -I/usr/include/python2.7 #-funroll loops -O2
LD = g++
LDFLAGS = -Wall -g -lm -std=c++11 -I/usr/include/python2.7 -lpython2.7
OBJS = utils.o prime-utils.o
	
default: gencerts.o $(OBJS)
	$(LD) -o $(EXE) $(OBJS) gencerts.o $(LDFLAGS) 
	
test: test.o $(OBJS)
	$(LD) -o test $(OBJS) test.o $(LDFLAGS) 

#install: default
#	install -c -m 0555 -o maxie -g bin $EXE ~/bin/

clean:
	rm -f *.o *.code $(EXE) test
	
clean-data: clean
	rm -f *.out

gencerts.o: prime-utils.h utils.h gencerts.cpp
	$(CC) -c $(CFLAGS) gencerts.cpp

prime-utils.o: prime-utils.h utils.h prime-utils.cpp
	$(CC) -c $(CFLAGS) prime-utils.cpp

utils.o: utils.h utils.cpp
	$(CC) -c $(CFLAGS) utils.cpp
	
test.o: utils.h utils.cpp prime-utils.h prime-utils.cpp test.cpp
	$(CC) -c $(CFLAGS) test.cpp

