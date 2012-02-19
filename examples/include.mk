CFLAGS+=-I../../thirdparty/samtools -I../../thirdparty/uthash/src -I../../inc -DUSE_BAM=1 -DHASH_FUNCTION=HASH_SFH
LDFLAGS+=-L../.. -L../../thirdparty/samtools -lkent -lbam -lm -lz
