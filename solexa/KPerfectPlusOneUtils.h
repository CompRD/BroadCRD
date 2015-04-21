#ifndef _K_PERFECT_PLUS_ONE_UTILS_H
#define _K_PERFECT_PLUS_ONE_UTILS_H

#include "String.h"

#include <cstdlib>
#include <unistd.h>

String getRealPath(String symlinkPath) {
	const char *filename = symlinkPath.c_str();
	char *buffer = NULL;
	int size = 64;
	int nchars;

	do {
		size *= 2;

		buffer = static_cast<char *>(realloc(buffer, size));
		nchars = readlink(filename, buffer, size);

		Assert(nchars != -1);
	} while (nchars == size);

	buffer[nchars] = '\0';
	String realPath(buffer);

	free(buffer);

	return realPath;
}

String getLookupTableName(String lookupPath) {
	return getRealPath(lookupPath).RevAfter("/").RevBefore(".fasta.lookuptable.lookup");
}

String getKTablePath(String UNIQUE_KMER_LOOKUP_PATH, String lookupPath, int K) {
	String ktablePath(UNIQUE_KMER_LOOKUP_PATH + "/" + getLookupTableName(lookupPath) + ".k" + ToString(K) + ".kp1table");

	return ktablePath;
}

String getKLookupPath(String HASHED_KMER_LOOKUP_PATH, String lookupPath, int K) {
	String klookupPath(HASHED_KMER_LOOKUP_PATH + "/" + getLookupTableName(lookupPath) + ".k" + ToString(K));

	return klookupPath;
}

#endif
