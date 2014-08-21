all:
	g++ -std=c++11 genomegen.cpp -o genomegen

clean:
	rm --force genomegen *.tab.* *.o *.lo test *.so*

cleantext:
	rm --force *.txt