///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "VecUtilities.h"
#include "util/RunCommand.h"

/**
 * class assembly_stats
 *
 * Container for assembly stats.
 */
class astats {

public:
  
  astats( const String *brief_file = 0 ) :
    name_( "" ),
    n_supers_( 0 ), n50_supers_( 0 ), tlen_supers_( 0 ),
    n_contigs_( 0 ), n50_contigs_( 0 ), tlen_contigs_( 0 ),
    covered_( -1 ), mutations_( -1 ), indels_( -1 )
  {
    if ( brief_file ) this->ParseBriefFile( *brief_file );
  }
  
  void ParseBriefFile( const String &brief_file ) {
    ifstream in( brief_file.c_str( ) );
    String key;
    String value;
    while ( in ) {
      in >> key >> value;
      if ( !in ) break;
      if ( key == "ASSEMBLY" ) name_ = value;
      else if ( key == "Nscaffolds" ) n_supers_ = value.Int( );
      else if ( key == "N50scaffolds" ) n50_supers_ = value.Int( );
      else if ( key == "TOTLENscaffolds" ) tlen_supers_ = value.Int( );
      else if ( key == "Ncontigs" ) n_contigs_ = value.Int( );
      else if ( key == "N50contigs" ) n50_contigs_ = value.Int( );
      else if ( key == "TOTLENcontigs" ) tlen_contigs_ = value.Int( );
      else if ( key == "COVERED" ) covered_ = value.Int( );
      else if ( key == "MUTATIONS" ) mutations_ = value.Int( );
      else if ( key == "INDELS" ) indels_ = value.Int( );
      else ForceAssert( 1 == 0 );
    }
  }
  
  String Name( ) const { return name_; }
  int NSupers( ) const { return n_supers_; }
  int N50Supers( ) const { return n50_supers_; }
  longlong TotLenSupers( ) const { return tlen_supers_; }
  int NContigs( ) const { return n_contigs_; }
  int N50Contigs( ) const { return n50_contigs_; }
  longlong TotLenContigs( ) const { return tlen_contigs_; }
  int Covered( ) const { return covered_; }
  int Mutations( ) const { return mutations_; }
  int Indels( ) const { return indels_; }

  void CompareWithOther( const astats *other, ostream &out ) const {
    vec< vec<String> > table;
    
    table.push_back( MkVec( String( "" ),
			    String( "new" ),
			    String( "base" ),
			    String( "diff (%)" ),
			    String( "" ) ) );

    this->CompareLength( true,
			 "n supers",
			 this->NSupers( ),
			 other->NSupers( ),
			 table );
    
    this->CompareLength( false,
			 "n50 supers",
			 this->N50Supers( ),
			 other->N50Supers( ),
			 table );
    
    this->CompareLength( false,
			 "tot length supers",
			 this->TotLenSupers( ),
			 other->TotLenSupers( ),
			 table );
    
    this->CompareLength( true,
			 "n contigs",
			this->NContigs( ),
			other->NContigs( ),
			table );

    this->CompareLength( false,
			 "n50 contigs",
			 this->N50Contigs( ),
			 other->N50Contigs( ),
			 table );
    
    this->CompareLength( false,
			 "tot length contigs",
			 this->TotLenContigs( ),
			 other->TotLenContigs( ),
			 table );
    
    this->CompareLength( false,
			 "covered",
			 this->Covered( ),
			 other->Covered( ),
			 table );
    
    this->CompareLength( true,
			 "mismatches",
			 this->Mutations( ),
			 other->Mutations( ),
			 table );

    this->CompareLength( true,
			 "indels",
			 this->Indels( ),
			 other->Indels( ),
			 table );

    this->CompareLength( true,
			 "tot errors",
			 this->Mutations( ) + this->Indels( ),
			 other->Mutations( ) + other->Indels( ),
			 table );

    out << "COMPARING\n"
	<< "  new:  " << this->Name( ) << "\n"
	<< "  base: " << other->Name( ) << "\n"
	<< "\n";
    PrintTabular( out, table, 3, "lrrrl" );
    out << endl;

  }


private:
  
  void CompareLength( const bool small_is_better,
		      const String &tag,
		      longlong newc,
		      longlong basec,
		      vec< vec<String> > &table ) const
  {
    String sdiff = "na";
    if ( basec != 0 ) {
      double pcdiff = 100. * SafeQuotient( newc - basec, basec );
      sdiff = ( newc > basec ? "+" : "" ) + ToString( pcdiff, 1 );
    }

    String descr = ".";
    if ( newc > basec ) descr = small_is_better ? "WORSE" : "better";
    else if ( newc < basec ) descr = small_is_better ? "better": "WORSE";
    
    table.push_back( MkVec( tag,
			    ToString( newc ),
			    ToString( basec ),
			    sdiff,
			    descr ) );
  }

  
private:

  String name_;            // name of assembly

  int n_supers_;           // super stats
  int n50_supers_;
  longlong tlen_supers_;

  int n_contigs_;          // contigs stats
  int n50_contigs_;
  longlong tlen_contigs_;

  longlong covered_;       // eval stats (optional)
  longlong mutations_;
  longlong indels_;
  
};

/**
 * CompareAssistedAssemblies
 *
 * Compare (assisted) assembly ASSEMBLY_NEW aginast ASSEMBLY_BASE, by
 * looking at the output of BuildAssistedReport.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY_NEW );
  CommandArgument_String( ASSEMBLY_BASE );
  EndCommandArguments;

  // File names.
  String new_file = ASSEMBLY_NEW + ".report.brief";
  String base_file = ASSEMBLY_BASE + ".report.brief";

  vec<String> needed;
  needed.push_back( new_file );
  needed.push_back( base_file );

  // Load assembly stats.
  astats new_stats( &new_file );
  astats base_stats( &base_file );

  // Compare.
  new_stats.CompareWithOther( &base_stats, cout );
  
}
