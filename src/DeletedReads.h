// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef DELETED_READS_H
#define DELETED_READS_H



/*
 * Cambridge, July 16, 2001
 *
 * DeletedReads.h
 *
 * Constants
 *  N_DELETION_CODES
 *  deletion_code
 *
 * Classes:
 *  del_datum
 *  deleted_reads
 *
 * Storage class for deleted reads.
 */
#include "String.h"
#include "Vec.h"



/*
 * Constants
 *
 * N_DELETION_CODES
 *  number of different deletion codes.
 *
 * If you change these codes, please also change CreateStandardOutputs.cc!
 *
 * deletion_code:
 *  code_config_noid: user defined (i.e.deleted from the configuration file;)
 *  code_trim_noid: not enough was left of it, after trimming;
 *    (Note: Reads which consist almost entirely (>= 85%) of one base
 *     and which have a quality score sum of <= 8000) are also included
 *     in this category.)
 *  code_match_vector: it matches E-coli or vector;
 *  code_mitochondrial: it or its mate appears to be mitochondrial;
 *  code_chimeric: it appears to be chimeric;
 *  code_multiple_instance: is 'multiple instance' of another read;
 *  code_other: read deleted for other reasons;
 *  code_unknown: read results unknown (to metainfo broker;)
 *  code_not_deleted: read has not been deleted.
 *  code_other_contaminant: read deleted by match to contaminants.fasta.
 *  code_same_name: delete all reads with the same name.
 *
 * WARNING: most but not all the reads in the *_noid categories have
 * been deleted before an id was assigned to them (see the note for
 * code_trim_noid).
 */
typedef char deletion_code;
const int N_DELETION_CODES = 11;

const char code_config_noid = 0;
const char code_trim_noid = 1;
const char code_match_vector = 2;
const char code_mitochondrial = 3;
const char code_chimeric = 4;
const char code_multiple_instance = 5;
const char code_other = 6;
const char code_unknown = 7;
const char code_not_deleted = 8;
const char code_other_contaminant = 9;
const char code_same_name = 10;



/*
 * Class del_datum
 *
 * Data for a deleted read. This is a triple of deletion code, deleted read
 * name, and deleted read id. Deletion code must always be given. Read name
 * may be "N/A", and read id may be -1, but one or the other must be given.
 * The idea is that sometimes the read name is available but the read id is
 * not, or vice versa.
 */ 
class del_datum {

public:

  del_datum(String read_name, int read_id, deletion_code del_code) :
    read_name_ ( read_name ), read_id_ ( read_id ), del_code_ ( del_code ) { }
  
  const String & GetName() const { return read_name_; }
  
  int GetReadId() const { return read_id_; }
  
  int GetDelCode() const { return (int)del_code_; }
  
  
private:
  
  String read_name_;  // Either deleted read name or "N/A".
  int read_id_;           // Either deleted read id or -1.
  int del_code_;          // Code for deletion.
};



/*
 * Class deleted_reads
 *
 * A simple class containing the names of deleted reads, with a
 * deletion reason. There are two files involved, containing
 * the deleted reads (with a deletion code,) and the unknown
 * reads (reads for which metainformations was not found.)
 *
 * Remark: must CreateFiles() before the first AddRead().
 */
class deleted_reads {

public:

  deleted_reads(const String &deleted_reads_file,
		const String &unknown_reads_file);
  
  // Creates new and empty "deleted reads" and "unknown" files. 
  void CreateFiles();

  void AddRead(const String& name, int id, deletion_code del_code) {
    if ( code_unknown == del_code )
      unk_stream_ << name << "\n";
    else
      del_stream_ << (int)del_code << "\t" << name << "\t" << id << "\n";
  }
  
  void AddRead(const String& name, deletion_code del_code) {
    AddRead( name, -1, del_code);
  }
  
  void AddRead(int id, deletion_code del_code) {
    AddRead( "N/A", id, del_code);
  }
  
  // Get deleted reads.
  vec<del_datum> GetDeletedReads();

  // Get unknown reads.
  vec<String> GetUnknownReads(); 

  // Parse deleted reads file and count how many "no ids" reads were deleted.
  int NoIdsCount( );

  // Find names of the deleted reads with no id. Return a sorted vector.
  void NoIdsNames( vec<String> &names );
  
  
private:
  
  String deleted_reads_file_;  // full path deleted reads file name;
  String unknown_reads_file_;  // full path unknown reads file name;
  ofstream del_stream_;        // ofstream fo deleted reads;
  ofstream unk_stream_;        // ofstream for unknown reads.
};



#endif
