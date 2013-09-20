## Run the program with `java TeX2Maxima`, which uses the .class file.
## However, you may have to add the classpath to your environment with
## something like the `export` statement below (or an equivalent
## `setenv` statement).

## This sets up java's "import" path
export CLASSPATH=.:snuggletex-1.2.2/bin/snuggletex-core-1.2.2.jar:snuggletex-1.2.2/bin/snuggletex-upconversion-1.2.2.jar:snuggletex-1.2.2/bin/saxon9-9.1.0.8.jar:snuggletex-1.2.2/bin/saxon9-dom-9.1.0.8.jar

all : TeX2Maxmima.class

TeX2Maxmima.class : TeX2Maxima.java
	@ javac TeX2Maxima.java

