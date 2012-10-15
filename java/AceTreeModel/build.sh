#!/bin/sh
# AceTree build script hack

# recompile
(cd src; ant)

# copy a bunch of text files which seem important.
# cp jar/org/rhwlab/snight/*.txt src/org/rhwlab/snight

# copy source files over
# cp -R src/* bin/

# jar command including old manifest.
# jar cfm AceTree_recomp.jar old_jar/META-INF/MANIFEST.MF -C bin .
jar cfm AceTree.jar src/META-INF/MANIFEST.MF -C src .

