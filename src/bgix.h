/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL), 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@live.co.uk> */

#ifndef __BGIX_H
#define __BGIX_H

#include <stdint.h>
#include "kstring.h"
#include "bgzf.h"

typedef int (*bi_fetch_f)(int l, const char *s, void *data);

struct __bi_index_t;
typedef struct __bi_index_t bi_index_t;

struct __bi_iter_t;
typedef struct __bi_iter_t *bi_iter_t;

typedef struct {
	BGZF *fp;
	bi_index_t *idx;
	char *fn, *fnidx;
} bgix_t;

typedef struct {
	int beg, end;
} bg_interval_t;

#ifdef __cplusplus
extern "C" {
#endif
	bgix_t *bi_open(const char *fn, const char *fnidx);
	int bi_lazy_index_load(bgix_t *t);
	void bi_close(bgix_t *t);
	bi_iter_t bi_query(bgix_t *t, int beg, int end);
	bi_iter_t ti_querys(bgix_t *t, const char *reg);
	const char *bi_read(bgix_t *t, bi_iter_t iter, int *len);

	/* Destroy the iterator */
	void bi_iter_destroy(bi_iter_t iter);

	/* Build the index for file <fn>. File <fn>.tbi will be generated
	 * and overwrite the file of the same name. Return -1 on failure. */
	int bi_index_build(const char *fn, const bi_conf_t *conf);

	/* Load the index from file <fn>.tbi. If <fn> is a URL and the index
	 * file is not in the working directory, <fn>.tbi will be
	 * downloaded. Return NULL on failure. */
	bi_index_t *bi_index_load(const char *fn);

	bi_index_t *bi_index_load_local(const char *fnidx);

	/* Destroy the index */
	void bi_index_destroy(bi_index_t *idx);

	/* Parse a region like: chr2, chr2:100, chr2:100-200. Return -1 on failure. */
	int bi_parse_region(const bi_index_t *idx, const char *str, int *begin, int *end);

	int bi_get_tid(const bi_index_t *idx, const char *name);

	/* Get the iterator pointing to the first record at the current file
	 * position. If the file is just openned, the iterator points to the
	 * first record in the file. */
	bi_iter_t bi_iter_first(void);

	/* Get the iterator pointing to the first record in region tid:beg-end */
	bi_iter_t bi_iter_query(const bi_index_t *idx, int beg, int end);

	/* Get the data line pointed by the iterator and iterate to the next record. */
	const char *bi_iter_read(BGZF *fp, bi_iter_t iter, int *len);

	const ti_conf_t *bi_get_conf(bi_index_t *idx);
	int bi_get_intv(const bi_conf_t *conf, int len, char *line, bi_interval_t *intv);

	/*******************
	 * Deprecated APIs *
	 *******************/

	/* The callback version for random access */
	int bi_fetch(BGZF *fp, const bi_index_t *idx, int tid, int beg, int end, void *data, bi_fetch_f func);

	/* Read one line. */
	int bi_readline(BGZF *fp, kstring_t *str);

#ifdef __cplusplus
}
#endif

#endif
