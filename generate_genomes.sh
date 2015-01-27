#!/bin/bash

# 10k
genomegen 10k 10 -s k --indel 0.002 --read-error 0.01 --read-garbage 0.05
zip 10k *.txt
rm -f *.txt

# 100k
genomegen 100k 100 -s k --indel 0.002 --read-error 0.01 --read-garbage 0.05
zip 100k *.txt
rm -f *.txt

# 1m
genomegen 1m 1 -s m --indel 0.002 --read-error 0.01 --read-garbage 0.05
zip 1m *.txt
rm -f *.txt

# 10m
genomegen 10m 10 -s m --indel 0.002 --read-error 0.01 --read-garbage 0.05
zip 10m *.txt
rm -f *.txt

# 100m
genomegen 100m 100 -s m --indel 0.002 --read-error 0.01 --read-garbage 0.05
zip 100m *.txt
rm -f *.txt

