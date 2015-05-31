#!/bin/sh

rm -r bin/bpcobar
JAVAFILES="src/bpcobar/BPCobar.java"
javac -Xlint -Xlint:-serial -source 6 -target 6 -bootclasspath lib/rt.jar -extdirs "" -d bin/ $JAVAFILES
jar cmf mainClass bpcobar.jar $JAVAFILES -C bin/ .

