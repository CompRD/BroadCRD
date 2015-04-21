///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
* \file Calculator.y
* \author ghall
* \date Oct 14, 2009
*
* \brief Simple infix calculator used by FeudalString (code adopted from the bison manual)
*/

// There is no target in the makefile to automatically parse yacc files, so
// if you edit this file, you will need to regenerate the code by hand.
//
//  $ bison Calculator.y -o Calculator.cc

// XXX: This is not reentrant
%{
#define YYSTYPE long
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
int yylex(void);
void yyerror(char const *);
const char *consume;
long result;
%}

/* Bison declarations.  */
%token      NUM
%left       '-' '+'
%left       '*' '/'
%left       NEG                     /* negation--unary minus */
%right      '^'                     /* exponentiation */

%% 
input:      /* empty */ 
            exp                     { result = $1; }
            ;

exp:        NUM                     { $$ = $1;          }
            | exp '+' exp           { $$ = $1 + $3;     } 
            | exp '-' exp           { $$ = $1 - $3;     }
            | exp '*' exp           { $$ = $1 * $3;     }
            | exp '/' exp           { $$ = $1 / $3;     }
            | '-' exp  %prec NEG    { $$ = -$2;         }
            | exp '^' exp           { $$ = pow($1, $3); }
            | '(' exp ')'           { $$ = $2;          }
            ;
%%

char _getchar()
{ char c = *consume; ++consume; return c; }

void _ungetc()
{ --consume; }

int yylex()
{ int c;
  char sz_chk[22];

  /* Skip white space. */ 
  while ((c = _getchar()) == ' ' || c == '\t')
    ;
  /* Process numbers. */
  if (isdigit(c))
  { _ungetc();
    sscanf(consume, "%ld", &yylval);
    sprintf(sz_chk, "%ld", yylval);
    consume += strlen(sz_chk);
    return NUM; }
  /* Return a single char. */
  return c;
}

void yyerror(char const *s)
{ fprintf(stderr, "%s\n", s); }

long evaluate(char const *s)
{ consume = s; yyparse(); return result; }

