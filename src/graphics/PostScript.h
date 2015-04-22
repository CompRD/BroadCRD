// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef POSTSCRIPT_H
#define POSTSCRIPT_H

#include <string>
#include <iostream>

#include "graphics/Color.h"

// This file contains tools for generating postscript.
//
// Namespace:
//  ns_psplot.

namespace ns_psplot {

  enum place { left, center, right };

  class freetext {
    
  public:
    
    freetext( const std::string &str, const color &col, const short font_size )
      : str_(str), col_(col), font_size_(font_size) { }
    
    std::string Str( ) const { return str_; }
    color Col( ) const { return col_; }
    short FontSize( ) const { return font_size_; }
    
    freetext( ) { }
    
  private:
    
    std::string str_;
    color col_;
    short font_size_;
    
  };
  
  
  inline std::string ps_string( const std::string &s ) { return " (" + s + ")"; }
  
  
  class text : public freetext {
    
  public:
    
    text( const std::string &str, const color &col, const short font_size, const place origin )
      : freetext( str, col, font_size ), origin_(origin) { }
    
    text( const freetext &f, const place origin ) : freetext(f), origin_(origin) { }
    
    text( ) { }
    
    friend std::ostream& operator<<( std::ostream& o, const text& t );
    
  private:
    
    place origin_;
    
  };


  class rlineto {
  public: 
    rlineto( float x, float y ) : x_(x), y_(y) { }
    friend std::ostream& operator<<( std::ostream& o, const rlineto& r )
      {
	return o << " " << r.x_ << " " << r.y_ << " rlineto";
      }
  private: 
    float x_, y_;
  };
  
  class lineto {
  public: 
    lineto( float x, float y ) : x_(x), y_(y) { }
    friend std::ostream& operator<<( std::ostream& o, const lineto& r )
      {
	return o << " " << r.x_ << " " << r.y_ << " lineto";
      }
  private: 
    float x_, y_;
  };
  
  class translate {
  public: 
    translate( float x, float y ) : x_(x), y_(y) { }
    friend std::ostream& operator<<( std::ostream& o, const translate& t )
      {  
	return o << " " << t.x_ << " " << t.y_ << " translate";
      }
  private: 
    float x_, y_;
  };
  
  class startpath {
  public: 
    startpath( float x, float y ) : x_(x), y_(y) { }
    friend std::ostream& operator<<( std::ostream& o, const startpath& s )
      {
	return o << " newpath " << s.x_ << " " << s.y_ << " moveto";
      }
  private: 
    float x_, y_;
  };
  
  class segment {
  public:  
    segment( float x1, float y1, float x2, float y2, const color &col )
      : x1_(x1), y1_(y1), x2_(x2), y2_(y2), col_(col) {}
    
    segment( ) { }
    
    float X1( ) const { return x1_; }
    
    float Y1( ) const { return y1_; }
    
    float X2( ) const { return x2_; }
    
    float Y2( ) const { return y2_; }
    
    float deltaX( ) const { return x2_ - x1_; }
    
    float deltaY( ) const { return y2_ - y1_; }
    
    color Color( ) const { return col_; }
    
    bool merge( const segment& next )
      {
	if ( x2_ != next.x1_ ||
	     y2_ != next.y1_ ||
	     col_ != next.col_ )
	  return false;
	else
	{
	  x2_ = next.x2_;
	  y2_ = next.y2_;
	  return true;
	}
      }
    
    friend std::ostream& operator<<( std::ostream& o, const segment& seg );
    
    friend bool operator==(const segment &seg1, const segment &seg2);
    
  private:
    float x1_, y1_, x2_, y2_;
    color col_;     
  };
  
}

#endif
