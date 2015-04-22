#ifndef _FASTBASEVECTOR_MAP_H
#define _FASTBASEVECTOR_MAP_H

#include "FastBasevector.h"
#include <ext/hash_map>
using __gnu_cxx::hash_map;

class fastbasevector_hash {
	public:
		size_t operator()(const fastbasevector &fb) const {
			const unsigned int _prime = 16777619;
			size_t _hash = 2166136261;

			for (unsigned int i = 0; i < fb.size(); ++i) {
				_hash ^= fb[i];
				_hash *= _prime;
			}

			return _hash;
		}
};

class fastbasevector_equal {
	public:
		bool operator()(const fastbasevector &fb1, const fastbasevector &fb2) const {
			return (fb1 == fb2);
		}
};

template<class T> class mapfastbasevector : public __gnu_cxx::hash_map<fastbasevector, T, fastbasevector_hash, fastbasevector_equal> {};

template<int K> class ptr_fast_base_t_hash {
    public:
        size_t operator()(fast_base_t* kmer) const {
            // djb2, from the internet.
            unsigned long hash = 5381;
            int n = 0;
            while (n < K) {
                hash = ((hash << 5) + hash) + *kmer; /* hash * 33 + c */
                kmer++;
                n++;
            }
            return hash;
        }
};

template<int K> class ptr_fast_base_t_equal {
    public:
        bool operator()(const fast_base_t* pfb1, const fast_base_t* pfb2) const {
            return strncmp(pfb1, pfb2, K) == 0;
        }
};

template<class T, int K> class mapptrfastbasevector : public __gnu_cxx::hash_map< fast_base_t *, T, ptr_fast_base_t_hash<K>, ptr_fast_base_t_equal<K> > {
    public:
        mapptrfastbasevector(unsigned int size) { this->resize(size); }
};

#endif
