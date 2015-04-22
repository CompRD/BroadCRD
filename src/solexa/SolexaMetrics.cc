/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "TokenizeString.h"
#include "solexa/SolexaMetrics.h"
#include "solexa/SolexaPipeline.h"

solexa_metric_db::solexa_metric_db( const filename_t& f, bool requireFileAlreadyExists)
{    if (!requireFileAlreadyExists && !IsRegularFile(f)) {
       return;
     }
     ReadMetrics( f );
}

Bool ExtractEntry( const String& x, int pos, String& answer )
{    if ( !x.Contains( "{", 0 ) || !x.Contains( "}", -1 ) ) return False;
     int commas = 0, lbrack = 0, rbrack = 0, last_comma = 0, n = x.size( );
     for ( int i = 1; i < n; i++ )
     {    if ( lbrack == rbrack && ( x[i] == ',' || i == n - 1 ) )
          {    ++commas;
               if ( commas == pos + 1 )
               {    answer = x.substr( last_comma + 1, i - last_comma - 1 );
                    return True;    }
               last_comma = i;    }
          if ( x[i] == '{' ) ++lbrack;
          if ( x[i] == '}' ) ++rbrack;    }
     return False;    }

String solexa_metric_db::Value( const String& key ) const
{    if ( key.Contains( "/" ) )
     {    String num = Value( key.Before( "/" ) ), denom = Value( key.After( "/" ) );
          return ToString( 100.0 * num.Double( ) / denom.Double( ), 1 ) + "%";    }
     maptype::const_iterator it = m_.find(key);
     if ( it == m_.end() )
     {    if ( key.Contains( "]", -1 ) )
          {    int pos1 = key.PosRev( "[" );
               if ( pos1 > 0 )
               {    String spos1 = key.substr( pos1+1, key.isize( ) - pos1 - 2 );
                    if ( spos1.IsInt( ) )
                    {    String key1 = key.substr( 0, pos1 );
                         maptype::const_iterator it1 = m_.find(key1);
                         if ( it1 != m_.end() )
                         {    String y1;
                              if ( ExtractEntry( it1->second, spos1.Int( ), y1 ) )
                                   return y1;    }
                         else if ( key1.Contains( "]", -1 ) )
                         {    int pos2 = key1.PosRev( "[" );
                              if ( pos2 > 0 )
                              {    String spos2 = key1.substr(
                                        pos2+1, key1.isize( ) - pos2 - 2 );
                                   if ( spos2.IsInt( ) )
                                   {    String key2 = key1.substr( 0, pos2 );
                                        maptype::const_iterator it2 = m_.find(key2);
                                        if ( it2 != m_.end() )
                                        {    String y1, y2;
                                             if ( ExtractEntry( it2->second,
                                                  spos2.Int( ), y1 ) )
                                             {    if ( ExtractEntry( y1,
                                                       spos1.Int( ), y2 ) )
                                                  {    return y2;    }    }    }
                                                      }    }    }    }    }    }
          FatalErr("Metric \"" << key << "\" is undefined.");    }
     return it->second;    }

bool solexa_metric_db::Defined( const String& key ) const
{    if ( key.Contains( "/" ) )
     {    String num = key.Before( "/" ), denom = key.After( "/" );
          return Defined(num) && Defined(denom);    }
     maptype::const_iterator it = m_.find(key);
     if ( it == m_.end() )
     {    if ( key.Contains( "]", -1 ) )
          {    int pos1 = key.PosRev( "[" );
               if ( pos1 > 0 )
               {    String spos1 = key.substr( pos1+1, key.isize( ) - pos1 - 2 );
                    if ( spos1.IsInt( ) )
                    {    String key1 = key.substr( 0, pos1 );
                         maptype::const_iterator it1 = m_.find(key1);
                         if ( it1 != m_.end() )
                         {    String y1;
                              if ( ExtractEntry( it1->second, spos1.Int( ), y1 ) )
                                   return true;    }
                         else if ( key1.Contains( "]", -1 ) )
                         {    int pos2 = key1.PosRev( "[" );
                              if ( pos2 > 0 )
                              {    String spos2 = key1.substr(
                                        pos2+1, key1.isize( ) - pos2 - 2 );
                                   if ( spos2.IsInt( ) )
                                   {    String key2 = key1.substr( 0, pos2 );
                                        maptype::const_iterator it2 = m_.find(key2);
                                        if ( it2 != m_.end() )
                                        {    String y1, y2;
                                             if ( ExtractEntry( it2->second,
                                                  spos2.Int( ), y1 ) )
                                             {    if ( ExtractEntry( y1,
                                                       spos1.Int( ), y2 ) )
                                                  {    return true;    }    }    }
                                                      }    }    }    }    }    }
          return false;    }
     return true;    }

bool solexa_metric_db::GetValue( const String& key, String& val ) const
{
  maptype::const_iterator it = m_.find(key);
  if (it==m_.end())
    return false;
  val = it->second;
  return true;
}

void solexa_metric_db::SetValue( const String& key, const String& val )
{
  m_[key] = val;
}

void solexa_metric_db::ReadMetrics( const filename_t& f ) {
  Ifstream( in, f );
  ReadMetrics( in, f );
}

void solexa_metric_db::ReadMetrics( ifstream& in, const filename_t& f, int maxLines )
{
  String line;
  int linesRead = 0;
  while(in &&  ( maxLines <= 0  ||  linesRead < maxLines )) {
    getline(in, line);
    if ( !in ) break;
    if ( line.Contains( "#", 0 ) ) continue;
    if ( !line.Contains( "=" ) ) {
      FatalErr("solexa_metric_db constructor: failed to parse file "
	       << f << "\n" << "problem line = \""
	       << line << "\"");
    }
    m_[line.Before( "=" )] = line.After( "=" );
    linesRead++;
  }
}

void solexa_metric_db::WriteMetrics( const filename_t& f, Bool append ) const
{
  ofstream out;
  if (append)
  { out.open( f.c_str( ), ios::app );
    if ( out.fail( ) ) FatalErr( "WriteMetrics: unable to open " << f << "." ); }
  else
  { out.open( f.c_str( ) );
    if ( out.fail( ) ) FatalErr( "WriteMetrics: unable to open " << f << "." ); }

  WriteMetrics(out);

  // Presumably, f = metric or <dir>/metric
  dirname_t targetDir;
  if (!f.Contains("/")) targetDir = ".";
  else targetDir = f.RevBefore("/");
  String theCommand = "touch " + targetDir;
  system(theCommand.c_str());

}

void solexa_metric_db::WriteMetricsMult( const filename_t& f, const vec< solexa_metric_db >& dbs ) {
  Ofstream( out, f );
  out << dbs.isize() << endl;
  For_( solexa_metric_db, db, dbs ) {
    out << db->m_.size() << endl;
    db->WriteMetrics( out );
  }
}


void solexa_metric_db::ReadMetricsMult( const filename_t& f, vec< solexa_metric_db >& dbs ) {
  Ifstream( in, f );
  String line;
  getline( in, line );
  ForceAssert( line.IsInt() );
  dbs.resize( line.Int() );
  for ( int i = 0; i < dbs.isize(); i++ ) {
    getline( in, line );
    ForceAssert( line.IsInt() );
    dbs[i].ReadMetrics( in, f, line.Int() );
  }
}


void solexa_metric_db::WriteMetrics( ostream& out ) const
{
  maptype::const_iterator it = m_.begin();
  for ( ; it!=m_.end(); ++it )
    out << it->first << "=" << it->second << "\n";
}
vec< vec<double> > solexa_metric_db::ValueVecVecDouble( const String& key ) const
{    String v = Value(key);
     ForceAssert( v.Contains( "{", 0 ) && v.Contains( "}", -1 ) );
     ForceAssert( !v.Contains( " " ) && !v.Contains( "\t" ) );
     v = v.After( "{" );
     v.resize( v.isize( ) - 1 );
     v.GlobalReplaceBy( "},{", "} {" );
     vec<String> vs;
     Tokenize( v, vs );
     vec< vec<double> > vv( vs.size( ) );
     for ( int i = 0; i < vs.isize( ); i++ )
         ParseDoubleSet(vs[i],vv[i],false);
     return vv;    }
