///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"

/**
 * RemoveContaminant
 *
 * Remove regions from contigs as specified by the given CONTAM file.
 * It updates contigs and scaffolds (superb), but it does not check
 * for scaffold integrity (linking-wise).
 *
 * CONTAM: full path name file (contig_id, contig_length, begin, end)
 * HEAD_IN: it loads <HEAD_IN>.{efasta,superb}
 * HEAD_OUT: it saves <HEAD_OUT>.{efasta,superb}
 * MIN_CONTIG_SIZE: minimum contig size to keep
 * JOIN_SEP_BY: simply join contigs separated by JOIN_SEP_BY or less
 */

Bool scompare( const superb& s1, const superb& s2 ){
  int l1 = s1.FullLength();
  int l2 = s2.FullLength();
  return (l1 > l2);
}

size_t scaffoldsTotLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].FullLength();
  return len;
}

size_t scaffoldsRedLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].ReducedLength();
  return len;
}

// restrict efasta given start position and length measured
// using the first option in the curly brackets.
void RestrictEfasta( efasta& source, size_t start, size_t len ){
  vec<size_t> enumer( source.size(), source.size() );
  size_t loc  = 0;
  for ( size_t i = 0; i < source.size(); i++ ){
    while( i < source.size() && source[i] != '{' ){
      enumer[i] = loc;
      loc++;
      i++;
     }
    while ( i < source.size() && source[i] != ',' ){    
      if ( source[i] != '{' ) 
	loc++;
       i++;    
    }
    if ( i < source.size() && source[i] == ',' ){
      size_t r = i;
      while( source[r] != '{' ){
	enumer[r] = loc;
	r--;
      }
      if ( source[r] == '{' )
	enumer[r] = loc;
    }
    while ( i < source.size() && source[i] != '}' ){
      enumer[i] = loc;
       i++;
    }
    if ( i < source.size() && source[i] == '}' )
      enumer[i] = loc;
  }    
 
  size_t startp = source.size();
  for ( int k = source.isize() -1; k >= 0; k-- )
    if ( enumer[k] == start )
      startp = k;

  size_t endp = source.size();
  for ( size_t k = 0; k < source.size(); k++ )
    if ( enumer[k] == start + len -1 )
      endp = k;
  endp++;

  ForceAssert( endp >= startp );
  ForceAssert( endp <= source.size() );
  
  source.erase( endp, source.size() - endp );
  source.erase( 0, startp );
  return;
}

// check integrity of scafolds and contigs data: contig size in superb == contig size in efasta,
//  each contig used once and only once
void check_integrity( const vec<superb>& scaffolds, 
		      const VecEFasta& efastas ){
  
  vec<int> tigs_used( efastas.size(), 0);
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas.size() );
      if ( efastas[tid].Length1() != s.Len(tpos) ){
	PRINT5( si, tpos, tid, s.Len(tpos), efastas[tid].Length1() );
	ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }
  vec<size_t> unused_tigs, overused_tigs;
  for ( size_t tid = 0; tid < tigs_used.size(); tid++ ){
    if ( tigs_used[tid] == 0 )
      unused_tigs.push_back( tid );
    else if ( tigs_used[tid] > 1 )
      overused_tigs.push_back(tid);
    
  }
  
  if ( unused_tigs.size() > 0 || overused_tigs.size() > 0 ){
    
    if ( unused_tigs.size() > 0 ){
      int max_un_size = efastas.at( unused_tigs[0] ).Length1();
      for ( size_t i = 0; i < unused_tigs.size(); i++ )
	if (  efastas.at( unused_tigs[i] ).Length1() > max_un_size )
	  max_un_size = efastas.at( unused_tigs[i] ).Length1();

      cout << "maximum size of unused contig is : " << max_un_size << endl;
    }

    PRINT2( unused_tigs.size(), overused_tigs.size() );
    ForceAssert( unused_tigs.size() == 0 && overused_tigs.size() == 0 );
  }
  
  return;
}



int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( CONTAM );
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Int_OrDefault_Doc( JOIN_SEP_BY, 0, 
      "Join contigs separated by JOIN_SEP_BY from each other in scaffolds");
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE, 1, 
      "Remove contigs that are smaller than MIN_CONTIG_SIZE");
  EndCommandArguments;

  // File names.
  String in_contam_file = CONTAM;
  String in_efasta_file = HEAD_IN + ".contigs.efasta";
  String in_superb_file = HEAD_IN + ".superb";
  String in_fasta_file  = HEAD_IN + ".contigs.fasta";
  String in_fastb_file  = HEAD_IN + ".contigs.fastb";
  
  String out_efasta_file = HEAD_OUT + ".contigs.efasta";
  String out_superb_file = HEAD_OUT + ".superb";
  String out_fastb_file  = HEAD_OUT + ".contigs.fastb";


  // reading contig information
  cout << Date( ) << " loading efasta file" << endl;
  if ( ! IsRegularFile( in_efasta_file) )
    FatalErr("input file " + in_efasta_file + " not found");
  VecEFasta efastas;
  LoadEfastaIntoStrings(in_efasta_file, efastas);
  cout << Date() << " converting efasta to basevectors" << endl;
 

  // Load contamination information
  cout << Date() << " loading contaminant information" << endl;
  if ( ! IsRegularFile( CONTAM ) )
    FatalErr("contamination file " + CONTAM + " not found");
  vec< vec<triple<int,int,int> > > v_contam( efastas.size() );
  int contamTotLen = 0;  // total contamination length
  ifstream in( CONTAM.c_str() );
  {
    String line;
    // reading data in the format:
    cout << Date() << "reading in contamination file" << endl;
    while ( getline(in,line) ){
      vec<String> tokens;
      Tokenize( line, tokens );
      ForceAssert( tokens.size() == 4);
      int tid    = tokens[0].Int();
      int tlen   = tokens[1].Int();
      int cbegin = tokens[2].Int();
      int cend   = tokens[3].Int();
      ForceAssertEq( tlen, efastas[tid].Length1() );
      ForceAssertGe( tid, 0 );
      ForceAssertLt( static_cast<size_t>(tid), efastas.size() );
      ForceAssertGe( cbegin, 0 );
      ForceAssertLe( cend, efastas[tid].Length1() );
      contamTotLen += cend - cbegin;
      v_contam.at(tid).push_back( triple<int,int,int>(tlen,cbegin, cend) );
    }
    in.close();
    cout << " Total contamination length = " << contamTotLen << endl;
  }


  // Lodading scaffolds
  cout << Date( ) << ": loading superb file" << endl;
  cout << " Original scaffold data" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( in_superb_file, scaffolds );
  
  size_t origScaffoldsTotLen = scaffoldsTotLen( scaffolds );
  size_t origScaffoldsRedLen = scaffoldsRedLen( scaffolds );


  // check initial basic integrity
  cout << Date() << " checking basic integrity" << endl;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssert( tid < efastas.size() );
      ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
    }
  }

  
  
  vec<superb> new_tscaffolds( efastas.size() ); 
  vec<Bool> modif_contigs( efastas.size(), False );
  vec< pair<int,int> > beg_gaps( efastas.size() );
  vec< pair<int,int> > end_gaps( efastas.size() );


  int contamCheckLen = 0;
  for ( size_t tid = 0; tid < v_contam.size(); tid++ ){  
    if ( v_contam.at(tid).size() == 0 ) 
      continue;
    else
      modif_contigs[tid] = True;
    cout << "\n  -------- \n\n"; PRINT(tid);
    int specContamLen = 0, resContamLen = 0;
    vec<Bool> treg( efastas[tid].Length1() +2, True );
    for ( size_t i = 0; i < v_contam[tid].size(); i++ ){
      ForceAssertEq( (int)efastas[tid].Length1(), v_contam[tid][i].first );
      specContamLen +=  v_contam[tid][i].third - v_contam[tid][i].second;
      for ( int j = v_contam[tid][i].second; j < v_contam[tid][i].third; j++ )
	treg[j+1] = False;
    }
    for ( size_t i = 1; i < treg.size() -1; i++ )
      if ( ! treg[i] ) {
	resContamLen++;
	contamCheckLen++;
      }
    ForceAssertEq( specContamLen, resContamLen );

    treg.front() = treg.back() = False; 

    vec<int> begs, ends;

    for ( size_t p = 1; p < treg.size(); p++ )
      if ( ! treg[p -1]  && treg[p] )
	begs.push_back( p - 1 );
      else if ( treg[p -1] && ! treg[p] )
	ends.push_back( p - 1 );
    
    ForceAssert( begs.size() == ends.size() );
    if ( begs.size() > 0 ){
      if ( begs.front() > 0 ){
	beg_gaps[tid].first  = begs[0];
	beg_gaps[tid].second = 1;
      }
      if ( ends.back() != efastas[tid].Length1() ){
	end_gaps[tid].first  = efastas[tid].Length1() - ends.back();
	end_gaps[tid].second = 1;
      }
    }else{
      // entire contig removed      
      int b = efastas[tid].Length1() / 2;
      int e = b + efastas[tid].Length1()%2;
      PRINT3( efastas[tid].Length1(), b, e );
      beg_gaps[tid] = pair<int,int>( b, 1 );
      end_gaps[tid] = pair<int,int>( e, 1 );
    }
    
    cout << "updating sequences" << endl;
    
    basevector tbases;
    efastas[tid].FlattenTo( tbases );
    efasta tefasta    = efastas[tid];
    new_tscaffolds[tid].SetNtigs( begs.size() );
    if ( begs.size() == 0 ){
      efastas[tid].resize(0);
    }
    for ( size_t i = 0; i < begs.size(); i++ ){
      basevector lbases( tbases, begs[i], ends[i] - begs[i] );
      efasta lefasta = tefasta;
      cout << "restricting efasta" << endl;
      RestrictEfasta( lefasta, begs[i], ends[i] - begs[i] );
      basevector tmpbases;
      lefasta.FlattenTo(tmpbases);
      ForceAssert( lbases == tmpbases );
      new_tscaffolds[tid].SetLen(i, lbases.size());
      cout << "updating efasta and base vecs" << endl;
      if ( i == 0 ){
	new_tscaffolds[tid].SetTig(i, tid);
	efastas[tid] = lefasta;
      }else{
	new_tscaffolds[tid].SetTig(i, efastas.size() );
	efastas.push_back( lefasta );
      }
      if ( i != begs.size() -1 ){
	new_tscaffolds[tid].SetGap(i, begs[i+1] - ends[i]);
	new_tscaffolds[tid].SetDev(i, 1);
      }
    }
  }
  if ( contamTotLen != contamCheckLen ){
    cout << " There seems to be a problem with contamination specification, are there contaminant overlaps?" << endl;
    ForceAssertEq( contamTotLen, contamCheckLen );
  }


  // update scaffolds
  size_t endGapSum = 0;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    superb & s = scaffolds[si];
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      if ( tid < modif_contigs.size() && modif_contigs[tid] ){
	cout << "\n ------------\n";
	new_tscaffolds[tid].PrintRange(cout, "insert:", 0, new_tscaffolds[tid].Ntigs() );
	size_t lenBefore = s.FullLength();
	size_t lenAfter = 0;
	if ( tpos == 0 ){
	  lenAfter += beg_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0 ){
	    lenAfter += end_gaps[tid].first;
	    lenAfter += s.Gap(tpos);
	  }
	}
	if ( tpos == s.Ntigs() -1 ) {
	  lenAfter += end_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0 ){
	    lenAfter += beg_gaps[tid].first;
	    lenAfter += s.Gap(tpos -1 );
	  }
	}
	endGapSum += lenAfter;
	s.PrintRange( cout, "before:", tpos -1, tpos +1 );
	s.ReplaceTigBySuper( tpos, new_tscaffolds[tid], 
			     beg_gaps[tid].first, beg_gaps[tid].second,
			     end_gaps[tid].first, end_gaps[tid].second);
	cout << "\nAfter" << endl;
	modif_contigs[tid] = False;
	s.PrintRange( cout, "after:", tpos -1, tpos + new_tscaffolds[tid].Ntigs() );
	lenAfter += s.FullLength();
	ForceAssertEq( lenBefore, lenAfter );
	tpos += new_tscaffolds[tid].Ntigs() -1;
      }
    }
  }  




  int joinedGapLen = 0;
  if ( JOIN_SEP_BY > 0 ){
    cout << "\n\n\n ------------------------\n" << endl;
    // joining very close contigs
    cout << Date() << " joining very close contigs" << endl;
    for ( size_t si = 0; si < scaffolds.size(); si++ ){
      superb & s = scaffolds[si];
      for ( int tpos = s.Ntigs() -1; tpos > 0; tpos-- ){
	size_t tid = s.Tig(tpos);
	if ( s.Gap(tpos -1 ) == JOIN_SEP_BY ){
	  joinedGapLen++;
	  size_t tidp = s.Tig(tpos -1);
	  efastas[tidp] = efastas[tidp] + efastas[tid];
	  efastas[tid].clear();
	  superb news;
	  news.SetNtigs(1);
	  news.SetTig( 0, tidp );
	  news.SetLen( 0, s.Len( tpos -1 ) + s.Len( tpos ) );
	  ForceAssertEq( efastas[tidp].Length1(), news.Len(0) );
	  news.SetTig( 0, tidp );
	  s.ReplaceTigsBySuper( tpos -1, tpos, news );
	}
      }
    }
    PRINT( joinedGapLen );
  }

  size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    newScaffoldsTotLen += scaffolds[is].FullLength();
    newScaffoldsRedLen += scaffolds[is].ReducedLength();
  }
  
  PRINT2( origScaffoldsTotLen, newScaffoldsTotLen + contamTotLen + joinedGapLen);
  PRINT2( origScaffoldsRedLen, newScaffoldsRedLen + contamTotLen );
  PRINT4( origScaffoldsTotLen, newScaffoldsTotLen, contamTotLen, joinedGapLen );
  PRINT3( origScaffoldsRedLen, newScaffoldsRedLen, contamTotLen );
  ForceAssertEq( origScaffoldsRedLen, newScaffoldsRedLen + contamTotLen );
  ForceAssertEq( origScaffoldsTotLen, newScaffoldsTotLen + endGapSum + joinedGapLen);
  
  {
    cout << Date() << " removing small contigs and renumbering" << endl;
    vec<int> offsets( efastas.size(), 0 );
    int offset = 0;
    for ( size_t tid = 0; tid < efastas.size(); tid++ ){
      if ( efastas[tid].Length1() < MIN_CONTIG_SIZE ){
	offsets[tid] = -1;
	offset++;
      }
      else{ offsets[tid] = offset; }
    }
    PRINT3( offsets.size(), efastas.size(), offset );
    
    ForceAssertEq( offsets.size(), efastas.size() );
    for ( int tid = 0; tid < offsets.isize(); tid++ ){
      if ( offsets[tid] > 0 )
	efastas[ tid - offsets[tid] ] = efastas[tid];
      
    }
    efastas.resize( efastas.size() - offset );
  
    cout << Date() << " updating scaffolds tig ids" << endl;
    for ( size_t si = 0; si < scaffolds.size(); si++ )
      for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
	int tid = scaffolds[si].Tig(tpos);
	if ( offsets[tid] >= 0 ){
	  int newtid = tid - offsets[tid];
	  ForceAssertGe( newtid, 0 );
	  scaffolds[si].SetTig( tpos, newtid );
	}
	else{
	  scaffolds[si].RemoveTigByPos( tpos );
	  tpos--;
	}
      }

  }


  cout << Date() << " removing empty scaffolds: " << endl;
  cout << "initial number of scaffolds = " << scaffolds.size() << endl;
  for ( int si = 0; si < scaffolds.isize(); si++ )
    if ( scaffolds.at(si).Ntigs() == 0 ){
      scaffolds.erase( scaffolds.begin() + si );
      si--;      
    }
  for ( int si = 0; si < scaffolds.isize(); si++ )
    ForceAssertGt( scaffolds.at(si).Ntigs(), 0 );
  cout << "final number of scaffolds = " << scaffolds.size() << endl;


  // Sorting scaffolds according to size and renumbering contigs according
  //   to sequential appearance in scaffolds
  cout << Date() << " sorting scaffolds" << endl;
  sort( scaffolds.begin(), scaffolds.end(), scompare);
 
  cout << Date() << " renumbering contigs for ordered scaffolds" << endl;
  vec<superb> oscaffolds = scaffolds;
  VecEFasta oefastas;
  int cTid = -1;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    for ( int tpos = 0; tpos < scaffolds[is].Ntigs(); tpos++ ){
      cTid++;
      int oTid = scaffolds[is].Tig(tpos);
      oscaffolds.at(is).SetTig( tpos, cTid );
      oefastas.push_back( efastas.at(oTid) );
    }
  }
  
  efastas.resize(0);

  // check integrity
  cout << Date() << " final integrity check" << endl;
  check_integrity( oscaffolds, oefastas );
  cout << Date() << " final integrity checked" << endl;

 
  // writing output
  cout << Date() << " writing superb" << endl;
  WriteSuperbs( out_superb_file, oscaffolds );

  cout << Date() << " writing efasta" << endl;
  Ofstream( efout, out_efasta_file );
  for ( size_t id = 0; id < oefastas.size(); id++ )
    oefastas[id].Print( efout, "contig_" + ToString(id) );

  vecbasevector obases( oefastas.size() );
  for ( size_t id = 0; id < oefastas.size(); id++ )
    oefastas[id].FlattenTo( obases[id] );
  
  obases.WriteAll( out_fastb_file );
  cout << Date() << " Done!" << endl;
}

