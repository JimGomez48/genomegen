all:
	g++ -std=c++11 genomegen.cpp -o genomegen

clean:
	rm --force genomegen *.tab.* *.o *.lo test *.so*

cleantext:
	rm --force *.txt

zip:
	zip genome ref_*.txt priv_*.txt ans_*.txt reads_*.txt
	rm ref_*.txt priv_*.txt ans_*.txt reads_*.txt

unzip:
	unzip genome.zip
	rm genome.zip
