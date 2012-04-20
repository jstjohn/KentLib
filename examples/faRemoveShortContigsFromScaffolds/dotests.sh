#!/bin/bash -x
./faRemoveShortContigsFromScaffolds -minGapLen=2 -minLen=200 scaffold_0.fa scaffold_0.trimmed.fa
./faRemoveShortContigsFromScaffolds -minGapLen=2 -minLen=200 test.fa test.out.fa
grep "AAGCATNCGGGAANNTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" *.fa
