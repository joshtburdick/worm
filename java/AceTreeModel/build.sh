#!/bin/sh
# AceTree build script hack

# recompile
(cd src; ant)

# copy a bunch of text files which seem important.
cp old_jar/org/rhwlab/snight/*.txt bin/org/rhwlab/snight

# jar command including old manifest.
jar cfm AceTree_recomp.jar old_jar/META-INF/MANIFEST.MF -C bin .

