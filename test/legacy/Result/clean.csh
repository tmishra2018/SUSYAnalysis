#!/bin/csh
find . -name "*.so" -exec rm -r {} \;
find . -name "*.d" -exec rm -r {} \;
find . -name "*.pcm" -exec rm -r {} \;
