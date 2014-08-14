all:
	g++ genomegen.cpp -o genomegen

clean:
	rm --force genomegen *.tab.* *.o *.lo test *.so*