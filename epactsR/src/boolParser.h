#ifndef __BOOL_PARSER_H
#define __BOOL_PARSER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <ctime>
#include <climits>

#define NUM_OP 5

// boolean parser accepting the following C-style expressions
// EVERY RULE has to written with whitespace
// ( expr1 && ( ! exprs2 ) ) || ( exprs3 ) ) && exprs4
class boolParser {
 public:
  std::vector<std::string> subs;
  std::vector<int> subl;
  std::vector<int> subm;
  //std::vector<bool> ans;
  std::vector<int> rules; // -1 : ( , -2 : ), -3 : &&, -4 : ||, -5 : !

  enum RULE { LBR = -1, RBR = -2, AND = -3, OR = -4, NOT = -5 };
  static const char* RULESTR[NUM_OP]; //= { "(", ")", "&&", "||", "!" };

  // example expression : ( ANNO=nonsynynmous:LDLR || ANNO=stopgain:LDLR ) && INDEL
  boolParser() {}

  // construct boolean parser
  boolParser(const char* s) {
    if ( s == NULL ) return;
    // srand(std::time(0)); // for debugging

    std::stringstream ss(s);
    std::string tok;
    for(int i=0; ss >> tok; ++i ) {
      if ( tok == "(" ) {
	rules.push_back(LBR);
      }
      else if ( tok == ")" ) {
	rules.push_back(RBR);
      }
      else if ( tok == "&&" ) {
	rules.push_back(AND);
      }
      else if ( tok == "||" ) {
	rules.push_back(OR);
      }
      else if ( tok == "!" ) {
	rules.push_back(NOT);
      }
      else {
	// check whether the operator is not separated by whitespace
	rules.push_back((int)subs.size());
	subs.push_back(tok);
	subl.push_back((int)tok.size());
	subm.push_back(0);
	//ans.push_back(false);

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
	o << subs[rules[i]];
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
	o << (subm[rules[i]] < 0 ? "TRUE" : "FALSE");
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

  bool parse(const char* s) {
    //print(std::cout,true);

    if ( rules.empty() ) return true;

    int ns = (int)subs.size();
    int t;
    for(t=0; t < ns; ++t) {
      subm[t] = 0;
    }
    char* p = (char*)s;
    int nm = 0;
    while( *p != '\0' ) {
      for(t=0; t < ns; ++t) {
	if ( subm[t] < 0 ) {
	  // do nothing
	}
	else if ( *p == subs[t][subm[t]] ) {
	  ++subm[t];
	  if ( subm[t] == subl[t] ) {
	    subm[t] = -1; // the substring was found
	    ++nm;
	    if ( nm == ns ) break;  // if all matches
	  }
	}
	else {
	  subm[t] = 0;
	}
      }
      ++p;
    }
    return eval();
  }

  bool parse(const char* s, int n) {
    //print(std::cout,true);

    if ( rules.empty() ) return true;

    int ns = (int)subs.size();
    int t;
    for(t=0; t < ns; ++t) {
      subm[t] = 0;
    }
    char* p = (char*)s;
    char* nch = p + n;
    int nm = 0;
    while( p < nch ) {
      for(t=0; t < ns; ++t) {
	if ( subm[t] < 0 ) {
	  // do nothing
	}
	else if ( *p == subs[t][subm[t]] ) {
	  ++subm[t];
	  if ( subm[t] == subl[t] ) {
	    subm[t] = -1; // the substring was found
	    ++nm;
	    if ( nm == ns ) break;
	  }
	}
	else {
	  subm[t] = 0;
	}
      }
      ++p;
    }
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
	int r = ( subm[rules[i]] < 0 ? 1 : 0 );
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
};

const char* boolParser::RULESTR[NUM_OP] = { "(", ")", "&&", "||", "!" };

#endif // __BOOL_PARSER_H
