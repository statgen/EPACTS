#ifndef __BOOL_EXPR_PARSER_H
#define __BOOL_EXPR_PARSER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <ctime>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <exception>

#define NUM_OP 11

class boolExprException : public std::exception {
public:
  virtual const char* what() const throw() {
    return "boolExprException was thrown";
  }
} beException;

// boolean expression parser accepting the following 
// EXAMPLE RULE :
// ( key1 && ( ! key2 ) ) || ( key3 >= -2 ) && ( key4 < 5 )
// EXAMPLE INPUT :
// key1;key2;key3=3.54;key4=2
class boolExprParser {
 public:
  // available operators
  enum BOOLOP { LBR = -1, RBR = -2, AND = -3, OR = -4, NOT = -5, GT = -6, GTE = -7, LT = -8, LTE = -9, EQ = -10, NE = -11 }; 
  static const char* RULESTR[NUM_OP];
  static double DBL_NAN;

  std::vector<int> rules;         // BOOLOP or index for keys
  std::vector<std::string> keys;  // substring keyword
  std::vector<int> lkeys;         // length of keyword
  std::vector<int> mkeys;         // match status of keyword
  std::vector<int> opkeys;        // operator of key
  std::vector<double> vkeys;      // value to be compared
  char* pEnd;

  // example expression : ( ANNO=nonsynynmous:LDLR || ANNO=stopgain:LDLR ) && INDEL
  // construct boolean parser
  boolExprParser(const char* s) {
    if ( s == NULL ) return;
    // srand(std::time(0)); // for debugging

    std::stringstream ss(s);
    std::string tok;
    for(int i=0; ss >> tok; ++i ) {
      // if known operator is observed, convert them as 'rules'
      if ( tok == "(" )       { rules.push_back(LBR); }
      else if ( tok == ")" )  { rules.push_back(RBR); }
      else if ( tok == "&&" ) { rules.push_back(AND); }
      else if ( tok == "||" ) { rules.push_back(OR);  }
      else if ( tok == "!" )  { rules.push_back(NOT); }
      else if ( tok == ">" )  { 
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = GT;  
      }
      else if ( tok == ">=" ) { 
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = GTE; 
      }
      else if ( tok == "<" )  { 
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = LT;  
      }
      else if ( tok == "<=" ) { 
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = LTE; 
      }
      else if ( tok == "==" ) {
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = EQ; 
      }
      else if ( tok == "!=" ) {
	if ( opkeys.back() < 0 ) throw beException;
	opkeys.back() = NE; 
      }
      else {
	if ( !opkeys.empty() && opkeys.back() < 0 && rules.back() >= 0 ) { 
	  double dtok = parse_numeric(tok);
	  if ( isnan(dtok) ) {
	    fprintf(stderr,"ERROR: Cannot parse %s\n",tok.c_str());
	    throw beException;
	  }
	  vkeys.back() = dtok;
	}
	else {
	  rules.push_back((int)keys.size());
	  keys.push_back(tok);
	  lkeys.push_back((int)tok.size());
	  mkeys.push_back(0);
	  opkeys.push_back(0);
	  vkeys.push_back(DBL_NAN);
	}

	// warn if operator is not delimited by whitespace
	for(int i=0; i < NUM_OP; ++i) {
	  if (tok.find(RULESTR[i]) != std::string::npos) {
	    std::cerr << "WARNING: Operator " << RULESTR[i] << " was found in the token " << tok << ". All operators must be delimited with whitespace" << std::endl;
	  }
	}
	//ans.push_back((rand()/(RAND_MAX+1.) < 0.5) ? true : false); // for debugging
      }
      //std::cerr << "#" << i << "\t" << rules[i] << "\t" << tok << std::endl;
    }
  }

  // print the current rule and values
  void print(std::ostream& o, bool runeval = true) {
    int n = (int)rules.size();
    o << "RULE:";
    for(int i=0; i < n; ++i) {
      o << "\t";
      if ( rules[i] < 0 ) {
	o << (const char*)RULESTR[0-(int)rules[i]-1];
      }
      else {
	o << keys[rules[i]];
      }
    }
    o << std::endl;

    o << "VALUE:";
    for(int i=0; i < n; ++i) {
      o << "\t";
      if ( rules[i] < 0 ) {
	o << RULESTR[0-rules[i]-1];
      }
      else {
	o << (mkeys[rules[i]] < 0 ? "TRUE" : "FALSE");
      }
    }
    if ( runeval ) o << "\t=\t" << eval();
    o << std::endl;
  }

  static void reduce(std::vector<int>& stack, int r) {
    while ( (! stack.empty()) && ( stack.back() != LBR ) ) {
      if ( stack.back() == NOT ) {
	stack.pop_back();
	r = 1-r;
      }
      else {
	int op = stack.back(); stack.pop_back(); // retrieve the operator
	int l = stack.back(); stack.pop_back();  // left operand
	if ( l < 0 ) {
	  std::cerr << "Error in the rule parsing : " << std::endl;
	  //print(std::cerr,false);
	  abort();
	}
	else {
	  if ( op == AND ) {
	    //stack.push_back((r && l) ? 1 : 0);
	    r = ((r && l) ? 1 : 0);
	  }
	  else if ( op == OR ) {
	    //stack.push_back((r || l) ? 1 : 0);
	    r = ((r || l) ? 1 : 0);
	    //std::cerr << "[ foo " << l << " " << r << " " << stack.back() << " ]" << std::endl;
	  }
	  else {
	    std::cerr << "Error in the rule parsing : " << std::endl;
	    //print(std::cerr,false);
	    abort();
	  }
	}
      }
    }
    stack.push_back(r);
  }

  // parse a rule string
  bool parse(const char* s, int n = -1) {
    // if no rule is specified, default is TRUE
    if ( rules.empty() ) return true;

    int ns = (int)keys.size();
    int t;
    for(t=0; t < ns; ++t) {
      mkeys[t] = 0;
    }
    char* p = (char*)s;
    char* nch = ( n < 0 ) ? NULL : p + n;
    int nm = 0;
    double strd;
    while( nch ? (p < nch) : (*p & 0xf0) ) { 
      for(t=0; t < ns; ++t) {
	if ( mkeys[t] < 0 ) {  // if already matches, do nothing
	}
	else if ( *p == keys[t][mkeys[t]] ) {  // if being matched, keep going on
	  ++mkeys[t];
	  if ( mkeys[t] == lkeys[t] ) { // substring match happened
	    if ( opkeys[t] < 0 ) { // if comparison operator is related
	      // check if the next character is '='
	      if ( p[1] == '=' ) {
		// parse the value
		strd = parse_numeric(p+2);
		switch(opkeys[t]) {
		case GT:
		  if ( strd > vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		case GTE:
		  if ( strd >= vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		case LT:
		  if ( strd < vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		case LTE:
		  if ( strd <= vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		case EQ:
		  if ( strd == vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		case NE:
		  if ( strd != vkeys[t] ) { mkeys[t] = -1; ++nm; }
		  else { mkeys[t] = 0; }
		  break;
		default:
		  fprintf(stderr,"Cannot recognize operator %d\n",opkeys[t]);
		  throw beException;
		}
	      }
	      else {
		mkeys[t] == 0; // declare mismatch
	      }
	    }
	    else {
	      mkeys[t] = -1; // the substring was found
	      ++nm;
	    }
	    if ( nm == ns ) break;  // if all matches
	  }
	}
	else {
	  mkeys[t] = 0;
	}
      }
      ++p;
    }
    pEnd = p;
    return eval();
  }

  // evaluate the boolean expression
  bool eval() {
    std::vector<int> stack;
    for(int i=0; i < (int)rules.size(); ++i) {
      //fprintf(stderr,"**** %d %d %s\n",i,rules[i],subs[i].c_str());

      switch(rules[i]) {
      case LBR:
      case AND:
      case OR:
      case NOT:
	stack.push_back(rules[i]); // simply insert the rule to stack before the next operand comes in
	break;
      case RBR: // search for the matching brace
	if ( stack.size() < 2 ) {
	  std::cerr << "Error in the rule parsing : empty between braces" << std::endl;
	  print(std::cerr,false);
	  abort();
	}
	else {
	  int v = stack.back(); stack.pop_back();   // value to return
	  int lbr = stack.back(); stack.pop_back(); // should match LBR
	  if ( ( lbr != LBR ) || ( v < 0 ) ) {
	    std::cerr << "Error in the rule parsing : " << std::endl;
	    print(std::cerr,false);
	    abort();
	  }
	  reduce(stack, v);
	}
	break;
      default: // resolve the rule
	int r = ( mkeys[rules[i]] < 0 ? 1 : 0 );
	reduce(stack, r);
      }
    }
    if ( stack.size() != 1 ) {
      std::cerr << "Error in the rule parsing : " << std::endl;
      print(std::cerr,false);
      abort();
    }
    else {
      return ( stack.back() > 0 ? true : false );
    }
  }

  static double parse_numeric(const std::string& str) {
    return parse_numeric(str.c_str(), (int)str.size());
  }

  static double parse_numeric(const char* s, int n) {
    char* pEnd;
    double d = strtod(s,&pEnd);
    if ( pEnd-s == n ) { return d; }
    else { return DBL_NAN; }
  }

  static double parse_numeric(const char* s) {
    char* pEnd;
    double d = strtod(s,&pEnd);
    if ( pEnd-s > 0 ) { return d; }
    else { return DBL_NAN; }
  }
};

const char* boolExprParser::RULESTR[NUM_OP] = { "(", ")", "&&", "||", "!", ">", ">=", "<", "<=", "==", "!=" };
double boolExprParser::DBL_NAN = sqrt(-1.);

#endif // __BOOL_EXPR_PARSER_H
