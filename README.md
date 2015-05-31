# bpcobar

Short and simple: given an element x in the cobar complex for BP_* BP, this program computes the coboundary dx.
For example: fire up the program and enter

    v2 [t1^2 t2 | t3] + v3
    
and you should see:

    -16 [ t1^7 ]
     + -4 [ t1t2^2 ]
     + 2 [ t3 ]
     + -56 v1 [ t1^6 ]
     + -1 v1 [ t2^2 ]
     + ... <61 more lines>

This program was written very quickly during collaborative computation sessions with Michael Andrews. It's not at all a polished program at this stage. If there are features or refinements that would be useful to see here, get in touch!


# How to use

Running the program requires only a single `jar` file, available at 

http://willperry.me/downloads/bpcobar-latest.jar

You'll want to run it from the command-line, e.g. by
    
    java -jar bpcobar-latest.jar



# Compiling

This repository includes a shell script `make` for compiling. Currently it's a bit specialized to my machine, in particular referring to a Java 6 runtime `rt.jar` in a subfolder, but it shouldn't be too hard to get it built if you're into Java. Feel free to contact me.


