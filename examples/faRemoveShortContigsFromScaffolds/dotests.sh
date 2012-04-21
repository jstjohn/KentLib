#!/bin/bash -x
echo "testing N splitting"
./faRemoveShortContigsFromScaffolds -minGapLen=2 -minLen=200 scaffold_0.fa scaffold_0.trimmed.fa
./faRemoveShortContigsFromScaffolds -minGapLen=2 -minLen=200 test.fa test.out.fa
./faRemoveShortContigsFromScaffolds -minGapLen=2 -minLen=200 test2.fa test2.out.fa
grep "NNTCNN" *.fa
echo "testing for smallest contig being too large"
./faRemoveShortContigsFromScaffolds -minGapLen=100 -minLen=200 scaffold_11934.fa scaffold_11934.trimmed.fa
./faSplitOn100Ns.pl  scaffold_11934.trimmed.fa >  scaffold_11934.trimmed.ctgs.fa
faSize scaffold_11934.trimmed.ctgs.fa
