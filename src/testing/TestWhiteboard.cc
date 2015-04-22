#include "MainTools.h"
#include "graphics/Whiteboard.h"

using namespace ns_whiteboard;

int main(int argc, char** argv )
{
  whiteboard aBoard;

  aBoard.Add( new point( make_pair(10,10) ) );
  aBoard.Add( new line( make_pair(20,20),
	    	       make_pair(30,30),
	    	       1.0, red ) );
  aBoard.Add( new text( make_pair(50,50),
			       "stuff", cyan ) );

  aBoard.Add( new arc( make_pair(50,50),
		       36.0,
		       0.0,
		       120.0,
		       1.0, green ) );

  ofstream psFile( "/tmp/doodad.eps" );

  ps_display aDisplay( psFile, 100, 100, 10 );

  aBoard.DisplayOn( &aDisplay );
}

