// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// alignment_caches are a way to avoid repeating the alignment of the
// same two reads over and over again.

// For a given (ordered) pair of ids, an RC value, and an offset, one
// can check the cache to see if the alignment has been done and
// retrieve the given alignment if so.  (Use the HasAlignment() method.)

// This is a template, so you can use alignment_cache<align>,
// alignment_cache<alignment>, alignment_cache<pack_align>, and so on.

template <class A>
class alignment_cache {
    
  private:

    // Class by which alignments are described.
    class alignment_key {
      public:
        alignment_key( const int id1, const int id2, 
                        const Bool rc2, const int offset )
            : id1_( id1 ), id2_( id2 ), rc2_( rc2 ), offset_( offset ) {}
        
        bool operator< ( const alignment_key &other ) const
        {
            if ( id1_ < other.id1_ ) return true;
            if ( id1_ > other.id1_ ) return false;
            if ( id2_ < other.id2_ ) return true;
            if ( id2_ > other.id2_ ) return false;
            if ( rc2_ < other.rc2_ ) return true;
            if ( rc2_ > other.rc2_ ) return false;
            if ( offset_ < other.offset_ ) return true;
            return false;
        }
        
      private:
        int id1_, id2_;
        Bool rc2_;
        int offset_;
    };
    
    // Class storing alignment and associated data.
    class alignment_data {
      public:
        alignment_data( const A &the_alignment, const float score )
            : the_alignment_( the_alignment ),               
              score_( score ) { }
        
        void Unpack( A &the_alignment, float &score ) const
        {
            the_alignment = the_alignment_;
            score = score_;
        }
        
      private:
        A the_alignment_;
        float score_;
    };

    // Class storing performance statistics of the cache.
    class alignment_cache_stats {
      public:
        alignment_cache_stats() 
            : checks_(0),
              hits_(0) {}
        
        void CacheHit()
        { 
            ++checks_;
            ++hits_;
        }
        
        void CacheMiss()
        {
            ++checks_;
        }

        void Print( ostream &out ) const
        {
            if ( checks_ == 0 )
                out << "No cache checks." << endl;
            else
                out << "Cache hit rate: " 
                    << PERCENT_RATIO( 3, hits_, checks_ ) << endl;
        }
        
      private:
        int checks_;
        int hits_;
    };
    
    // This is the main data store for the cache.
    map<alignment_key,alignment_data> aligns_;
    
    // Stats for the cache.  We use a pointer to avoid const conflicts.
    alignment_cache_stats *stats_;
    
  public:
    alignment_cache()
        : stats_( new alignment_cache_stats ) {}
    
    ~alignment_cache()
    {
        delete stats_;
    }


    // AddAlignment()

    // Add an alignment to the cache.

    void AddAlignment( const int id1, const int id2, 
                       const Bool rc2, const int offset,
                       const A &the_alignment, const float score )
    {
        alignment_key key( id1, id2, rc2, offset );
        alignment_data data( the_alignment, score );
        aligns_.insert( make_pair( key, data ) );
    }
    

    // HasAlignment()

    // Check if the cache has an alignment for the given id1, id2,
    // rc2, and offset.  If there is such an alignment, fill out the
    // score and the_alignment parameters with the cached values and
    // return true.  If there is no such alignment, return false.

    bool HasAlignment( const int id1, const int id2,
                       const Bool rc2, const int offset,
                       A &the_alignment, float &score ) const
    {
        typename map<alignment_key,alignment_data>::const_iterator align_iter;
        
        alignment_key target( id1, id2, rc2, offset );
        align_iter = aligns_.find( target );
        
        if ( align_iter == aligns_.end() )
        {
            stats_->CacheMiss();
            return false;
        }

        stats_->CacheHit();
        align_iter->second.Unpack( the_alignment, score );
        return true;
    }


    // PrintStats()

    // Print performance statistics of the cache.

    void PrintStats( ostream &out ) const
    {
        out << "Aligns cached: " << aligns_.size() << endl;
        stats_->Print( out );
    }
};
