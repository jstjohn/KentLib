#include <string.h>
#include <assert.h>
#include "kstring.h"
#include "bam.h"

static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
	if (n != b->core.n_cigar) {
		int o = b->core.l_qname + b->core.n_cigar * 4;
		if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data) {
			b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
			kroundup32(b->m_data);
			b->data = (uint8_t*)realloc(b->data, b->m_data);
		}
		memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
		memcpy(b->data + b->core.l_qname, cigar, n * 4);
		b->data_len += (n - b->core.n_cigar) * 4;
		b->core.n_cigar = n;
	} else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

#define write_cigar(_c, _n, _m, _v) do { \
		if (_n == _m) { \
			_m = _m? _m<<1 : 4; \
			_c = (uint32_t*)realloc(_c, _m * 4); \
		} \
		_c[_n++] = (_v); \
	} while (0)

static void unpad_seq(bam1_t *b, kstring_t *s)
{
	int k, j, i;
	uint32_t *cigar = bam1_cigar(b);
	uint8_t *seq = bam1_seq(b);
	ks_resize(s, b->core.l_qseq);
	for (k = 0, s->l = 0, j = 0; k < b->core.n_cigar; ++k) {
		int op, ol;
		op = bam_cigar_op(cigar[k]);
		ol = bam_cigar_oplen(cigar[k]);
		assert(op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CSOFT_CLIP);
		if (op == BAM_CMATCH) {
			for (i = 0; i < ol; ++i) s->s[s->l++] = bam1_seqi(seq, j);
			++j;
		} else if (op == BAM_CSOFT_CLIP) {
			j += ol;
		} else {
			for (i = 0; i < ol; ++i) s->s[s->l++] = 0;
		}
	}
}

int bam_pad2unpad(bamFile in, bamFile out)
{
	bam_header_t *h;
	bam1_t *b;
	kstring_t r, q;
	uint32_t *cigar2 = 0;
	int n2 = 0, m2 = 0, *posmap = 0;

	h = bam_header_read(in);
	bam_header_write(out, h);
	b = bam_init1();
	r.l = r.m = q.l = q.m = 0; r.s = q.s = 0;
	while (bam_read1(in, b) >= 0) {
		uint32_t *cigar = bam1_cigar(b);
		n2 = 0;
		if (b->core.pos == 0 && b->core.tid >= 0 && strcmp(bam1_qname(b), h->target_name[b->core.tid]) == 0) {
			int i, k;
			unpad_seq(b, &r);
			write_cigar(cigar2, n2, m2, bam_cigar_gen(b->core.l_qseq, BAM_CMATCH));
			replace_cigar(b, n2, cigar2);
			posmap = realloc(posmap, r.m * sizeof(int));
			for (i = k = 0; i < r.l; ++i) {
				posmap[i] = k; // note that a read should NOT start at a padding
				if (r.s[i]) ++k;
			}
		} else {
			int i, k, op;
			unpad_seq(b, &q);
			if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) write_cigar(cigar2, n2, m2, cigar[0]);
			for (i = 0, k = b->core.pos; i < q.l; ++i, ++k)
				q.s[i] = q.s[i]? (r.s[k]? BAM_CMATCH : BAM_CINS) : (r.s[k]? BAM_CDEL : BAM_CPAD);
			for (i = k = 1, op = q.s[0]; i < q.l; ++i) {
				if (op != q.s[i]) {
					write_cigar(cigar2, n2, m2, bam_cigar_gen(k, op));
					op = q.s[i]; k = 1;
				} else ++k;
			}
			write_cigar(cigar2, n2, m2, bam_cigar_gen(k, op));
			if (bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CSOFT_CLIP) write_cigar(cigar2, n2, m2, cigar[b->core.n_cigar-1]);
			for (i = 2; i < n2; ++i)
				if (bam_cigar_op(cigar2[i]) == BAM_CMATCH && bam_cigar_op(cigar2[i-1]) == BAM_CPAD && bam_cigar_op(cigar2[i-2]) == BAM_CMATCH)
					cigar2[i] += cigar2[i-2], cigar2[i-2] = cigar2[i-1] = 0;
			for (i = k = 0; i < n2; ++i)
				if (cigar2[i]) cigar2[k++] = cigar2[i];
			n2 = k;
			replace_cigar(b, n2, cigar2);
			b->core.pos = posmap[b->core.pos];
		}
		bam_write1(out, b);
	}
	free(r.s); free(q.s); free(posmap);
	bam_destroy1(b);
	bam_header_destroy(h);
	return 0;
}

int main_pad2unpad(int argc, char *argv[])
{
	bamFile in, out;
	if (argc == 1) {
		fprintf(stderr, "Usage: samtools depad <in.bam>\n");
		return 1;
	}
	in = strcmp(argv[1], "-")? bam_open(argv[1], "r") : bam_dopen(fileno(stdin), "r");
	out = bam_dopen(fileno(stdout), "w");
	bam_pad2unpad(in, out);
	bam_close(in); bam_close(out);
	return 0;
}
