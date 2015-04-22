#ifndef _FAST_BASEVECTOR_H
#define _FAST_BASEVECTOR_H

#include "Basevector.h"
#include "Vec.h"
#include "String.h"

typedef char fast_base_t;
const fast_base_t FAST_BASE_A = 0x0d;
const fast_base_t FAST_BASE_T = 0x0c;
const fast_base_t FAST_BASE_G = 0x0b;
const fast_base_t FAST_BASE_C = 0x0a;
const fast_base_t FAST_BASE_UNKNOWN = 0x00;

class fastbasevector {
	public:
		fastbasevector() : _block_size(0), _block(NULL) {};
		fastbasevector(unsigned int n, int extra=0) : _block_size(0), _block(NULL) { allocate(n+extra); };
		fastbasevector(String &kmer) : _block_size(0), _block(NULL) {
			initializeFromString(kmer);
		}

		fastbasevector(char const*ckmer) : _block_size(0), _block(NULL) {
			String skmer(ckmer);

			initializeFromString(skmer);
		}

		fastbasevector(const fastbasevector &fb) : _block_size(0), _block(NULL) {
			allocate(fb.size());

			memcpy(_block, fb.GetBlockPtr(), _block_size);
		}

		~fastbasevector() { deallocate(); }

		static fast_base_t GetComplementaryBase(fast_base_t base) {
			return base^1;
		}

		static char as_base(fast_base_t base) {
			switch (base) {
				case FAST_BASE_A: return 'A'; break;
				case FAST_BASE_C: return 'C'; break;
				case FAST_BASE_G: return 'G'; break;
				case FAST_BASE_T: return 'T'; break;
				default:          return '?'; break;
			};
		}

		fast_base_t operator[](unsigned int i) const {
			return _block[i];
		}

		friend bool operator==(const fastbasevector &fb1, const fastbasevector &fb2) {
			if (fb1._block_size != fb2._block_size) { return 0; }

			return (!memcmp(fb1._block, fb2._block, fb1._block_size));
		}

		friend bool operator!=(const fastbasevector &fb1, const fastbasevector &fb2) {
			return !(fb1 == fb2);
		}

		fastbasevector& operator=(const fastbasevector &fb) {
			if (&fb != this) {
				allocate(fb.size());

				memcpy(_block, fb.GetBlockPtr(), _block_size);
			}

			return *this;
		}

		fastbasevector& operator=(basevector &b) {
			allocate(b.size());

			for (unsigned int i = 0; i < _block_size; ++i) {
				base_t base = b[i];

				switch (base) {
					case BASE_A: _block[i] = FAST_BASE_A; break;
					case BASE_C: _block[i] = FAST_BASE_C; break;
					case BASE_G: _block[i] = FAST_BASE_G; break;
					case BASE_T: _block[i] = FAST_BASE_T; break;
				};
			}

			return *this;
		}

		friend bool operator<(const fastbasevector &x, const fastbasevector &y) {
			for (unsigned int i = 0; i < x.size(); ++i) {
				if (i >= y.size()) { return false; }
				if (x[i] < y[i])   { return true;  }
				if (x[i] > y[i])   { return false; }
			}
			return x.size() < y.size();
		}

		unsigned int size(void) const { return _block_size; }
		int isize(void) const { return static_cast<int>(_block_size); }

		fast_base_t *GetBlockPtr(unsigned int i = 0) const { return &_block[i]; }

		void SetToSubOf(const fastbasevector &orig_fbv, unsigned int start_pos, int len, int extra=0) {
			allocate(len);

			memcpy(_block, orig_fbv.GetBlockPtr(start_pos), len);
		}

		void SetToSubOf(const basevector &orig_bv, unsigned int start_pos, int len, int extra=0) {
			allocate(len);

			for (int i = 0; i < len; ++i) {
				base_t base = orig_bv[start_pos + i];

				switch (base) {
					case BASE_A: _block[i] = FAST_BASE_A; break;
					case BASE_C: _block[i] = FAST_BASE_C; break;
					case BASE_G: _block[i] = FAST_BASE_G; break;
					case BASE_T: _block[i] = FAST_BASE_T; break;
					default:     _block[i] = FAST_BASE_UNKNOWN; break;
				};
			}
		}

		void ReverseComplement(void) {
			fast_base_t temp1, temp2;

			for (unsigned int i = 0; i < _block_size/2; ++i) {
				temp1 = _block[i];
				temp2 = _block[_block_size - i - 1];
				_block[i] = (temp2^1);
				_block[_block_size - i - 1] = (temp1^1);
			}

			if (_block_size % 2 != 0) { _block[(_block_size/2)] ^= 1; }
		}

		String ToString(void) const {
			String s;
			s.resize(_block_size);

			for (unsigned int i = 0; i < _block_size; ++i) {
				s[i] = fastbasevector::as_base(_block[i]);
			}

			return s;
		}

	private:
		void allocate(unsigned int size) {
			_block_size = size;
			_block = reinterpret_cast<fast_base_t *>(realloc(_block, _block_size*sizeof(fast_base_t)));

			Assert(_block != NULL);
		}

		void deallocate(void) {
			if (_block != NULL) {
				free(_block);
				_block = NULL;
				_block_size = 0;
			}
		}

		void initializeFromString(String &kmer) {
			allocate(kmer.size());

			for (unsigned int i = 0; i < kmer.size(); ++i) {
				switch (kmer[i]) {
					case 'A':
					case 'a': _block[i] = FAST_BASE_A; break;

					case 'C':
					case 'c': _block[i] = FAST_BASE_C; break;

					case 'G':
					case 'g': _block[i] = FAST_BASE_G; break;

					case 'T':
					case 't': _block[i] = FAST_BASE_T; break;

					default:  _block[i] = FAST_BASE_UNKNOWN; break;
				};
			}
		}

		unsigned int _block_size;
		fast_base_t *_block;
};

class vecfastbasevector {
	public:
		vecfastbasevector(String filename) {
			vecbasevector vbv(filename);

			unsigned int size = vbv.size();
			fbvs.resize(size);

			for (unsigned int i = 0; i < size; ++i) {
				fbvs[i] = vbv[i];
			}
		}

		unsigned int size(void) { return fbvs.size(); }

		fastbasevector &operator[](unsigned int i) { return fbvs[i]; }

	private:
		vec<fastbasevector> fbvs;
};

#endif
