////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
//////////////////////////////////////////////////////////////////////////// 

#ifndef GENOME_INFO_CACHE_H
#define GENOME_INFO_CACHE_H

#include "Basevector.h"
#include "Intvector.h"
#include "String.h"
#include "Vec.h"

template <typename T>
class GenomeInfoCache
{
    private:
    vec<vecbasevector::size_type> block_indices;
    vec<vecbasevector::size_type> intra_block_indices;
    vecbasevector::size_type current_block_index;
    vec<vecbasevector::size_type> first_block_chr;
    vec<vecbasevector::size_type> chr_per_block;
    String cache_dir;
    vec<T> current_block;
    static const size_t MAX_BLOCK_SIZE = 100000000;
    bool erase_at_destruction;
    bool read_only;

    public:
    GenomeInfoCache(const vecbasevector& genome, const String& dir,
                    bool is_temporary = true, bool read_only_cache = false,
                    bool in_memory = false);
    T* getInfo(vecbasevector::size_type chr);
    T& operator[](vecbasevector::size_type chr);
    void setReadOnly(bool ro_state);
    size_t size();
    ~GenomeInfoCache();

    private:
    void readCache(vecbasevector::size_type cache_index);
    void updateCache();
    void constructBlock(vecbasevector::size_type cache_index);
    String blockFilename(vecbasevector::size_type cache_index);
};

template <typename T>
GenomeInfoCache<T>::GenomeInfoCache(const vecbasevector& genome,
                                    const String& dir, bool is_temporary,
                                    bool read_only_cache, bool in_memory) :
    erase_at_destruction(is_temporary), read_only(read_only_cache)
{
    block_indices.resize(genome.size(), 0);
    intra_block_indices.resize(genome.size(), 0);
    current_block_index = 0;

    size_t block_size_counter = 0;
    size_t num_blocks = 0;
    size_t ib_ind = 0;
    for (size_t chr = 0; chr < genome.size(); chr++)
    {
        if (ib_ind == 0)
        {
            first_block_chr.push_back(chr);
            chr_per_block.push_back(0);
        }
        
        block_indices[chr] = num_blocks;
        block_size_counter += genome[chr].size();
        intra_block_indices[chr] = ib_ind;
        ib_ind++;
        chr_per_block[num_blocks]++;

        if (!in_memory && block_size_counter > MAX_BLOCK_SIZE)
        {
            num_blocks++;
            block_size_counter = 0;
            ib_ind = 0;
        }
    }

    if (is_temporary)
    {
        if (num_blocks > 1)
        {
            cache_dir = dir + "/cov_XXXXXX";
            char* res = mkdtemp(const_cast<char*>(cache_dir.c_str()));
            if (res == NULL)
            {
                std::cerr << "failed to create directory " << cache_dir
                    << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
            std::cerr << "created cache directory " << cache_dir << std::endl;
        }
        else // no need for using the disk at all
        {
            constructBlock(0);
        }
    }
    else
    {
        cache_dir = dir;
#ifdef DEBUG
        std::cerr << "GENOME_INFO_CACHE DEBUG: set directory to " << dir
            << std::endl;
#endif // DEBUG
    }

    return;
}

template <typename T>
T* GenomeInfoCache<T>::getInfo(vecbasevector::size_type chr)
{
    // load the cached block, if necessary
    if (current_block.size() == 0 || current_block_index != block_indices[chr])
    {
        readCache(block_indices[chr]);
    }
#ifdef DEBUG
    std::cerr << "GENOME_INFO_CACHE DEBUG: mapping " << chr << " to "
         << intra_block_indices[chr] << std::endl;
#endif // DEBUG
    return &(current_block[intra_block_indices[chr]]);
}

template <typename T>
T& GenomeInfoCache<T>::operator[](vecbasevector::size_type chr)
{
    return *getInfo(chr);
}

template <typename T>
void GenomeInfoCache<T>::setReadOnly(bool ro_state)
{
    if (!read_only && ro_state)
    {
        updateCache();
    }
    read_only = ro_state;
    return;
}

template <typename T>
size_t GenomeInfoCache<T>::size()
{
    return block_indices.size();
}

template <typename T>
void GenomeInfoCache<T>::updateCache()
{
    String cache_file = blockFilename(current_block_index);
#ifdef DEBUG
    std::cerr << "GENOME_INFO_CACHE DEBUG: writing " << cache_file;
#endif // DEBUG
    BinaryWriter::writeFile(cache_file.c_str(), current_block);
#ifdef DEBUG
    std::cerr << " done" << std::endl;
#endif // DEBUG
}

template <typename T>
void GenomeInfoCache<T>::readCache(vecbasevector::size_type cache_index)
{
    // save the block if it exists
    if (current_block.size() != 0 && cache_index != current_block_index
        && !read_only)
    {
        updateCache();
    }

    String new_block_cache_file = blockFilename(cache_index);

    if (IsRegularFile(new_block_cache_file))
    {
#ifdef DEBUG
        std::cerr << "GENOME_INFO_CACHE DEBUG: reading "
            << new_block_cache_file;
#endif // DEBUG
        current_block.clear();
        BinaryReader::readFile(new_block_cache_file.c_str(), &current_block);
#ifdef DEBUG
        std::cerr << " done" << std::endl;
#endif // DEBUG
    }
    else // block not previously requested, construct from scratch
    {
        constructBlock(cache_index);
    }

    current_block_index = cache_index;

    return;
}

template <typename T>
void GenomeInfoCache<T>::constructBlock(vecbasevector::size_type cache_index)
{
#ifdef DEBUG
    std::cerr << "GENOME_INFO_CACHE DEBUG: constructing block " << cache_index;
#endif // DEBUG
    current_block.clear();
    current_block.resize(chr_per_block[cache_index]);
#ifdef DEBUG
    std::cerr << " done" << std::endl;
#endif // DEBUG
    return;
}

template <typename T>
String GenomeInfoCache<T>::blockFilename(vecbasevector::size_type cache_index)
{
    return cache_dir + "/block_" + ToString(cache_index);
}

template <typename T>
GenomeInfoCache<T>::~GenomeInfoCache()
{
    if (erase_at_destruction)
    {
        for (size_t bi = 0; bi <= block_indices[block_indices.size() - 1]; bi++)
        {
            String block_cache_file = blockFilename(bi);
            if (IsRegularFile(block_cache_file))
            {
#ifdef DEBUG
                std::cerr << "GENOME_INFO_CACHE DEBUG: deleting "
                     << block_cache_file;
#endif // DEBUG
                if (unlink(block_cache_file.c_str()) != 0)
                {
                    std::cerr << "GenomeInfoCache destructor failed to delete "
                         << block_cache_file << std::endl;
                }
#ifdef DEBUG
                std::cerr << " done" << std::endl;
#endif // DEBUG
            }
        }

        if (IsDirectory(cache_dir))
        {
#ifdef DEBUG
            std::cerr << "GENOME_INFO_CACHE DEBUG: deleting " << cache_dir;
#endif // DEBUG
            if (rmdir(cache_dir.c_str()) != 0)
            {
                std::cerr << "failed to remove " << cache_dir << std::endl;
            }
#ifdef DEBUG
            std::cerr << " done" << std::endl;
#endif // DEBUG
        }
    }
    else if (!read_only)
    {
        updateCache();
    }
    return;
}

#endif // GENOME_INFO_CACHE_H
