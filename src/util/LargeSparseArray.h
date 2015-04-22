/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are
 * reserved.

 * This software is supplied without any warranty or guaranteed support
 * whatsoever. Neither the Broad Institute nor MIT can be responsible for its
 * use, misuse, or functionality.

 * LargeSparseArray.h
 * Efficiently implements a human genome sized array.
 */

#ifndef LARGE_SPARSE_ARRAY_H_
#define LARGE_SPARSE_ARRAY_H_

#include "system/Assert.h"
#include "Basevector.h"
#include "CommonSemanticTypes.h"
#include "Vec.h"
#include "system/file/FileReader.h"
#include <fstream>
#include <list>


clock_t cheap_clock()
{
	static clock_t ticks = 0;
	ticks += 1; 
	return ticks;
}

// It is assumed that RECORD can be initialized with memset(&record, 0, sizeof(RECORD),
// and that RECORD has no ctor nor dtor
template <class RECORD>
class LargeSparseArray
{

public:
    struct Header
    {
        unsigned long magic_number;
        unsigned long total_array_size;
        unsigned long block_size;
        unsigned long number_of_blocks;
    };

private:
    struct TOCRecord
    {
        unsigned long block_id;
        unsigned long file_offset;

        TOCRecord() { block_id = 0; file_offset = 0; }
    };
    
    class CacheRecord 
    {
    public:
      unsigned long block_id;
      clock_t age;
      CacheRecord(unsigned long block_id):
        block_id(block_id),
        age(cheap_clock())
        {}
    };

    struct BlockWithCacheIndex
    {
      static const unsigned long NULL_CACHE_ID = ULONG_MAX;
      // Points to RECORD[block_size] elements 
      RECORD *ptr;
      // Index of this block in the cache, or NULL_CACHE_ID if this block is not in cache.
      unsigned long cache_id;
    };
    
    static const clock_t CLOCK_MAX = LONG_MAX;

    Header         header;
    // table_of_contents.size() == header.total_array_size
    // blocks.size() == header.total_array_size
    vec<TOCRecord> table_of_contents;
    vec<BlockWithCacheIndex>   blocks;

    unsigned long cache_size;
    // cache.size() == cache_size
    vec<CacheRecord>  cache;

    unsigned long file_size;
    std::fstream file;

    bool read_only;
    
    RECORD null_record;

    // Ensures that blocks[block_id] points to an allocated block
    void get_a_page(unsigned long block_id)
    {
      Assert(blocks[block_id].ptr == NULL);
      if (cache.size() < cache_size) {
        blocks[block_id].ptr = (RECORD *)new char[sizeof(RECORD) * header.block_size];
        blocks[block_id].cache_id = cache.size();
        cache.push_back(CacheRecord(block_id));
        return;
      }
      // we've got to drop the oldest page.
      unsigned long min_index = cache.size();
      clock_t       min_value = CLOCK_MAX;
      unsigned int  i;
      for (i = 0; i < cache.size(); i++)
      {
        //                    printf(">> %d %d %d %d\n", i, min_index, min_value, cache_ages[i]); fflush(stdout);

        if (cache[i].age < min_value)
        {
          //                       printf(">>! %d %d %d %d\n", i, min_index, min_value, cache_ages[i]); fflush(stdout);
          min_index = i;
          min_value = cache[i].age;
        }
      }

      //                printf("freeing block %lu (0x%x)\n", min_index, blocks[min_index]); fflush(stdout);

      if (!read_only) { sync_block(cache[min_index].block_id); file.flush(); }
      unsigned long old_block_id = cache[min_index].block_id;
      RECORD *theBlock = blocks[old_block_id].ptr;
      blocks[old_block_id].ptr = NULL;
      blocks[old_block_id].cache_id = BlockWithCacheIndex::NULL_CACHE_ID;
      blocks[block_id].ptr = theBlock;
      blocks[block_id].cache_id = min_index;
      cache[min_index].block_id = block_id;
      cache[min_index].age = cheap_clock();
    }
    
    bool block_allocated(unsigned long block_id) {
      if (blocks[block_id].ptr != NULL) {
        return true;
      }
      return (table_of_contents[block_id].file_offset != 0);
    }
    
    void touch_page(unsigned long block_id) 
    {
      Assert(blocks[block_id].ptr != NULL);
      cache[blocks[block_id].cache_id].age = cheap_clock();
    }

    bool fetch(unsigned long block_id, bool allocate=false)
    {
        bool verbose=false;
        if (verbose)
        {
	        printf("(fetch %lu %lu ",
	                    block_id,
	                    table_of_contents[block_id].file_offset);
	        fflush(stdout);
        }

        if (table_of_contents[block_id].file_offset == 0)
        {
            // this block isn't allocated yet.
            if (verbose) { printf("not_allocated)\n"); fflush(stdout); }

            if (allocate == false) { return false; }

            get_a_page(block_id);
            memset(blocks[block_id].ptr, 0x00, sizeof(RECORD)*header.block_size);
            table_of_contents[block_id].block_id    = block_id;
            table_of_contents[block_id].file_offset = file_size;
            file_size += header.block_size * sizeof(RECORD);
            return true;
        }

        else if (blocks[block_id].ptr == NULL)
        {
            // this block exists, but isn't loaded.
            if (verbose) { printf("not_loaded)\n"); fflush(stdout); }

            get_a_page(block_id);

            memset(blocks[block_id].ptr, 0x00, sizeof(RECORD)*header.block_size);
            file.seekg(table_of_contents[block_id].file_offset, std::ios::beg);
            file.read((char*)blocks[block_id].ptr, sizeof(RECORD)*header.block_size);

            return true;
        }

        // this block exists and is loaded, so do nothing.
        else
        {
            if (verbose) { printf("hit)\n"); fflush(stdout); }
            touch_page(block_id);
            return true;
        }
    }

    void sync_header()
    {
        ForceAssert(this->read_only == false);

        file.seekp( 0x00, std::ios::beg );
        file.write((char const*)&header, sizeof(Header));
        for (unsigned long block_id = 0; block_id < blocks.size(); block_id++)
        {
            file.write((char const*)&table_of_contents[block_id], sizeof(TOCRecord));
        }
    }

    void sync_block(unsigned long block_id)
    {
        ForceAssert(this->read_only == false);

//        printf("syncing block %lu\n", block_id); fflush(stdout);

        if ((table_of_contents[block_id].file_offset != 0) && (blocks[block_id].ptr != NULL))
        {
            file.seekp(table_of_contents[block_id].file_offset, std::ios::beg);
            file.write((char const*)blocks[block_id].ptr, sizeof(RECORD)*header.block_size);
        }

        //fflush(file);
    }

private:
  // Disallow copy ctor and operator=
  LargeSparseArray(const LargeSparseArray&);
  LargeSparseArray &operator=(const LargeSparseArray&);
    
public:

    LargeSparseArray(unsigned long cache_size=50000):
      cache_size(cache_size),
      read_only(true)
      {
        memset(&null_record, 0, sizeof(RECORD));
      }

    ~LargeSparseArray(void)
    {
        for (unsigned long block_id = 0; block_id < blocks.size(); block_id++)
        {
            delete [] blocks[block_id].ptr;
        }
        close();
    }
    
    unsigned int blocks_loaded() const {
      return cache.size();
    }

    unsigned long get_block_size() const 
    {
        return header.block_size;
    }

    unsigned long get_total_array_size() const 
    {
        return header.total_array_size;
    }
    

    void set(unsigned long offset, RECORD r)
    {
        ForceAssert(this->read_only == false);

        unsigned long block_id = offset / header.block_size;
        if ((memcmp(&r, &null_record, sizeof(RECORD)) == 0) && !block_allocated(block_id)) {
          // No need to do anything if writing nulll value into non-instantiated block
          return;
        }
        fetch(block_id, true);
        unsigned long block_offset = offset - (block_id*header.block_size);
        //printf("(set %lu %lu %lu 0x%lx)\n", offset, block_id, block_offset, blocks[block_id]); fflush(stdout);
        blocks[block_id].ptr[block_offset] = r;
    }

    bool get(unsigned long offset, RECORD& r, bool allocate=false)
    {
        if (allocate)
        {
            ForceAssert(this->read_only==false);
        }

        unsigned long block_id = offset / header.block_size;
        bool status = fetch(block_id, allocate);
        unsigned long block_offset = offset - (block_id*header.block_size);

        if (status == false) { memset(&r, 0x00, sizeof(RECORD)); return false; }
        else                 { r = blocks[block_id].ptr[block_offset];  return true;  }
    }

    RECORD* get(unsigned long offset, bool allocate=false)
    {
        if (allocate)
        {
            ForceAssert(this->read_only==false);
        }

        unsigned long block_id = offset / header.block_size;
        bool status = fetch(block_id, allocate);
        unsigned long block_offset = offset - (block_id*header.block_size);
        if (status == false) { return NULL; }
        else                 { return &(blocks[block_id].ptr[block_offset]); }
    }

    static void readHeader(Header *header, const String &file_name) {
      FileReader fr(file_name.c_str());
      fr.read(header,sizeof(Header));
    }
    
    void open(String file_name, bool read_only=true)
    {
        this->read_only = read_only;

        struct stat s;
        stat(file_name.c_str(), &s);
        file_size = s.st_size;

        file.open(file_name.c_str(),
                    read_only ?
                            std::ios::in|std::ios::binary :
                            std::ios::in|std::ios::binary|std::ios::out);
        ForceAssert(file);

        file.seekg(0, std::ios::beg);
        file.read((char*)&header, sizeof(Header));

        table_of_contents.resize(header.number_of_blocks);
        blocks.resize(header.number_of_blocks);
        Assert(cache.empty());
        cache.reserve(cache_size);

        for (unsigned long block_id = 0; block_id < blocks.size(); block_id++)
        {
            file.read((char*)&table_of_contents[block_id], sizeof(TOCRecord));
        }
    }

    void create(String file_name, unsigned long block_size, unsigned long size)
    {
        this->read_only=false;

        file.open(file_name.c_str(),
                   std::ios::in|std::ios::out|std::ios::binary|std::ios::trunc);
        ForceAssert(file);

        header.magic_number     = 0xDEADBEEF;
        header.total_array_size = size;
        header.block_size       = block_size;
        header.number_of_blocks = (size / block_size) + 1;

        table_of_contents.resize(header.number_of_blocks);
        blocks.resize(header.number_of_blocks);
        Assert(cache.empty());
        cache.reserve(cache_size);

        file_size = sizeof(Header) + (sizeof(TOCRecord) * table_of_contents.size());

        file.seekp(0, std::ios::beg);
        file.write((char const*)&header, sizeof(Header));
        for (unsigned long block_id = 0; block_id < blocks.size(); block_id++)
        {
            file.write((char const*)&table_of_contents[block_id], sizeof(TOCRecord));
        }

    }

    void load_all(bool allocate=false)
    {
        this->cache_size = table_of_contents.size();
        cache.reserve(this->cache_size);
        for (unsigned int i = 0; i < table_of_contents.size(); i++)
        {
            fetch(i, allocate);
        }
    }

    void sync()
    {
        ForceAssert(this->read_only == false);


        printf("syncing %lu blocks out of %lu total\n", cache.size(), blocks.size()); 
        fflush(stdout);

        sync_header();
        for (unsigned long i = 0; i < cache.size(); ++i)
        {
            sync_block(cache[i].block_id);
        }
        file.flush();
    }

    void close()
    {
      if (file.is_open()) {
        file.close();
      }
    }

    unsigned long size() { return header.total_array_size; }

    void report()
    {
        printf("LargeSparseArrayReport : %lu %lu %lu %lu %lu\n",
                   header.total_array_size,
                   header.block_size,
                   header.number_of_blocks,
                   header.number_loaded,
                   header.number_loaded * header.block_size * sizeof(RECORD));
        fflush(stdout);
    }

};

size_t genomic_to_linear(int contig, int offset, vecbasevector& ref)
{
	size_t ans = 0;
	for (size_t i = 0; i < ref.size(); i++)
	{
		if (i != static_cast<size_t>(contig)) { ans += ref[i].size(); }
		else { ans += offset; break; }
	}
	return ans;
}

void linear_to_genomic(size_t linear, int& contig, int& offset, vecbasevector& ref)
{
    size_t original_linear = linear;
    for (contig = 0; static_cast<size_t>(contig) < ref.size(); contig++)
    {
        if (linear < ref[contig].size()) { offset = linear; return; }
        else { linear -= ref[contig].size(); }
    }
    printf("ERROR linear_to_genomic: original linear = %lu\n", original_linear);
    ForceAssert(0);
}

size_t reference_length(vecbasevector& ref)
{
    return genomic_to_linear(ref.size()-1, ref[ref.size()-1].size()-1, ref);
}

#endif /* LARGE_SPARSE_ARRAY_H_ */
