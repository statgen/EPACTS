/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Parameters.h"
#include "Constant.h"
//#include "MathConstant.h"
#include "Error.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctype.h>
#include <stdarg.h>

int Parameter::nameCol = 30;
int Parameter::statusCol = 15;

Parameter::Parameter(char c, const char * desc, void * v)
{
  ch = (char) tolower(c);
  description = new char [strlen(desc) + 1];
  strcpy(description, desc);
  var = v;
  warnings = NULL;
}

bool Parameter::Read(int , char ** argv, int argn)
{
  int   p = 0;
  char  c = (char) tolower(argv[argn][p]);

  if ((c == '-') || (c == '/'))
    {
      p++;
      c = (char) tolower(argv[argn][p]);
    }

  if (c == ch)
    {
      Translate(&(argv[argn][++p]));
      return true;
    }
  return false;
}

bool Parameter::TranslateExtras(const char * , const char *)
{
  return false;
}

void Parameter::warning(const char * format, ...)
{
  va_list ap;
  va_start(ap, format);
  char buf[65535];
  vsprintf(buf, format, ap);
  va_end(ap);

  if (warnings == NULL) {
    ::warning(buf);
  }
  else {
    (*warnings) += buf;
  }
}

void IntParameter::Translate(const char * value)
{
  //fprintf(stderr,"--- %s\n",value);
  *(int *) var = atoi(value);
}

bool IntParameter::TranslateExtras(const char * value, const char * extras)
{
  //fprintf(stderr,"--- %s %s\n",value, extras);

  if (value[0] != 0 || !CheckInteger(extras))
    return false;

  Translate(extras);

  return true;
}

void IntParameter::Status()
{
  fprintf(stderr, "%*s : %*d (-%c9999)\n", nameCol, description,
	  statusCol, *(int *) var, ch);
}

void SwitchParameter::Translate(const char * value)
{
  notice("SwitchParameter::Translate(%s)",value);
  switch (*value)
    {
    case '+' :
      *(bool *) var = true;
      break;
    case '-' :
      *(bool *) var = false;
      break;
    case 0 :
      *(bool *) var = ! * (bool *) var;
      break;
    default :
      warning("Command line parameter -%c%s: the option '%c' has no meaning\n", ch, value, value[0]);
    }
}

void SwitchParameter::Status()
{
  fprintf(stderr, "%*s : %*s (-%c[+|-])\n", nameCol, description,
	  statusCol, *(bool *) var == false ? "OFF" : "ON", ch);
}

DoubleParameter::DoubleParameter(char c, const char * desc, double & v)
  : Parameter(c, desc, &v)
{
  precision = 2;
}

void DoubleParameter::Translate(const char * value)
{
  if (value[0])
    *(double *) var = atof(value);
  else
    *(double *) var = _NAN_;
}

bool DoubleParameter::TranslateExtras(const char * value, const char * extras)
{
  if (value[0] != 0 || !CheckDouble(extras))
    return false;

  Translate(extras);

  return true;
}

void DoubleParameter::Status()
{
  double absolute_value = fabs(* (double *) var);

  if (*(double *) var == _NAN_)
    fprintf(stderr, "%*s : %*s (-%c99.999)\n", nameCol, description,
	    statusCol, "NAN", ch);
  else if (absolute_value >= 0.00095)
    fprintf(stderr, "%*s : % *.*f (-%c99.999)\n", nameCol, description,
	    statusCol, precision, * (double *) var, ch);
  else if (absolute_value <= 1e-15)
    fprintf(stderr, "%*s : % *.0f (-%c99.999)\n", nameCol, description,
	    statusCol, * (double *) var, ch);
  else
    fprintf(stderr, "%*s : %*.0e (-%c99.999)\n", nameCol, description,
	    statusCol, *(double *) var, ch);
}

void StringParameter::Translate(const char * value)
{
  std::string * s = (std::string *) var;

  *s = value;
}

bool StringParameter::TranslateExtras(const char * value, const char * extras)
{
  if ((value[0] != 0) || ((!required) && (extras[0] == '-')))
    return false;

  std::string * s = (std::string *) var;

  *s = extras;

  return true;
}

void StringParameter::Status()
{
  fprintf(stderr, "%*s : %*s (-%cname)\n", nameCol, description,
	  statusCol, ((std::string *) var)->c_str(), ch);
}

void ListParameter::Status()
{
  OptionList * l;

  for (l = options; l->ch != 0; l++)
    if (l->code == *((int *)var))
      break;

  fprintf(stderr, "%*s : %*s (-%c[%s])\n", nameCol, description,
	  statusCol, l->description, ch, key.c_str());
}

void ListParameter::Translate(const char * value)
{
  OptionList * l;

  for (l = options; l->ch != 0; l++)
    if (tolower(l->ch) == tolower(value[0]))
      break;

  if (l->ch == 0 && tolower(value[0]) != 0)
    warning("Command line parameter -%c%s: the option '%c' has no meaning\n",
	    ch, value, value[0], key.c_str());

  *((int*) var) = l->code;
}

ListParameter::ListParameter(char c, const char * desc, int & v, OptionList * opt)
  : Parameter(c, desc, &v)
{
  options = opt;

  for (OptionList * l = options; l->ch != 0; l++)
    {
      key += l->ch;
      key += '|';
    }

  key.resize((int)key.size() - 1);
}

SetParameter::SetParameter(char c, const char * desc, int & v, OptionList * opt)
  : Parameter(c, desc, &v)
{
  options = opt;

  for (OptionList * l = options; l->ch != 0; l++)
    {
      key += l->ch;
      key += '|';
    }
  key.resize(key.size() - 1);
}

void SetParameter::Status()
{
  bool first = 0;
  int  temp = * (int *) var;

  for (OptionList * l = options; l->ch != 0; l++)
    if ((l->code & temp) || (l->code == *(int *) var))
      {
	if (!first)
	  fprintf(stderr, "%*s : %*s (-%c{%s})\n", nameCol, description,
		  statusCol, l->description, ch, key.c_str());
	else
	  fprintf(stderr, "%*s & %*s\n", nameCol, "",
		  statusCol, l->description);
	first = true;
	temp &= ~l->code;
      }
}

void SetParameter::Translate(const char * value)
{
  *(int*)var = 0;

  for (const char * chr = value; *chr != 0; chr++)
    {
      int valid = false;

      for (OptionList * l = options; l->ch != 0; l++)
	if (tolower(l->ch) == tolower(*chr))
	  {
	    *((int*) var) |= l->code;
	    valid = true;
	  }

      if (!valid)
	warning("Command line parameter -%c%s: the option '%c' has no meaning\n",
		ch, value, *chr);
    }
}

LongParameters::LongParameters(const char * desc, LongParameterList * lst)
  : Parameter('-', desc, NULL)
{
  list = lst;

  index.clear();
  legacyIndex.clear();
  group_len = 0;

  LongParameterList * ptr = list + 1;

  while (ptr->description != NULL)
    {
      if (ptr->type == LP_LEGACY_PARAMETERS)
	break;

      if (ptr->value != NULL)
	index[ptr->description] = ptr;
      else {
	int tmp = strlen(ptr->description);
	if (tmp > group_len) group_len = tmp;
      }

      ptr++;
    }

  while (ptr->description != NULL)
    {
      if (ptr->value != NULL)
	legacyIndex[ptr->description] = ptr;

      ptr++;
    }

  precision = 2;
}

/*
void LongParameters::ExplainAmbiguity(const char * cstr)
{
  std::string value(cstr);

  int p = value.FastFindChar(':');
  std::string stem = p == -1 ? value : value.substr(0,p);
  std::string matches;

  for(StringMapIterator i = index.begin(); i != index.end(); ++i) {
    if ( i->first == stem ) {
      if (matches.size() + i->first.size() > 50) {
	matches += " ...";
	break;
      }
      matches += " --";
      matches += i->first;
    }
  }

  warning("Ambiguous --%s matches%s\n",
	  value.c_str(), matches.c_str());
}
*/

// Translate long parameters, now allows only exact matches
void LongParameters::Translate(const char * cstr)
{
  std::string value(cstr);
  size_t p = value.find(':');
  std::string query = (p == value.npos) ? value : value.substr(0,p);
  StringMapIterator i = index.find(query);

  //fprintf(stderr,"*** cstr = %s\n",cstr);

  LongParameterList * ptr;
  if ( i != index.end() ) {
    ptr = (LongParameterList*) i->second;
    //fprintf(stderr,"*** cstr = %s\n",cstr);
  }
  else {

    i = legacyIndex.find(query);
    if (i == legacyIndex.end())
    {
      //fprintf(stderr,"Command line parameter --%s is undefined\n", value.c_str());
      error("Command line parameter --%s is undefined\n", value.c_str());
      return;
    }
    ptr = (LongParameterList *) i->second;
    ptr->touched = true;
  }

  if (ptr->type == LP_BOOL_PARAMETER) {
    if (p == value.npos) {
      * (bool *) ptr->value ^= true;
    }
    else
	* (bool *) ptr->value = ((value.substr(p+1) == "ON") ? true : false);
    
    // In exclusive groups, only one option may be selected
    if (ptr->exclusive)
      {
	for (int i = -1; ptr[i].exclusive; i--) *(bool *)ptr[i].value = false;
	for (int i =  1; ptr[i].exclusive; i++) *(bool *)ptr[i].value = false;
      }
  }
  else if (ptr->type == LP_INT_PARAMETER) {
    //fprintf(stderr,"*** %s\n",value.c_str());
    if (p == value.npos )
      * (int *) ptr->value = * (int *) ptr->value ? 0 : 1;
    else
      *(int *) ptr->value = value.substr(p + 1) == "ON" ?
	1 : atoi(value.substr(p + 1).c_str());
  }
  else if (ptr->type == LP_DOUBLE_PARAMETER)
    {
      if (p != value.npos)
	* (double *) ptr->value = atof(value.substr(p + 1).c_str());
    }
  else if (ptr->type == LP_STRING_PARAMETER)
    {
      if (p != value.npos)
	* (std::string *) ptr->value = value.substr(p + 1);
    }
}

    /*
  String value(cstr);

  int p = value.FastFindChar(':');
  int option = p == -1 ? index.FindStem(value) : index.FindStem(value.Left(p));

  if (option == -2)
    {
      ExplainAmbiguity(cstr);
      return;
    }

  LongParameterList * ptr;

  if (option >= 0)
    ptr = (LongParameterList *) index.Object(option);
  else
    {
      int alternate = p == -1 ? legacyIndex.FindFirstStem(value) :
	legacyIndex.FindFirstStem(value.Left(p));

      if (alternate < 0)
	{
	  warning("Command line parameter --%s is undefined\n", (const char *) value);
	  return;
	}

      ptr = (LongParameterList *) legacyIndex.Object(alternate);
      ptr->touched = true;
    }

  if (ptr->type == LP_BOOL_PARAMETER)
    {
      if (p == -1)
	* (bool *) ptr->value ^= true;
      else
	*(bool *) ptr->value = value.SubStr(p + 1).SlowCompare("ON") == 0;

      // In exclusive groups, only one option may be selected
      if (ptr->exclusive)
	{
	  for (int i = -1; ptr[i].exclusive; i--) *(bool *)ptr[i].value = false;
	  for (int i =  1; ptr[i].exclusive; i++) *(bool *)ptr[i].value = false;
	}
    }
  else if (ptr->type == LP_INT_PARAMETER)
    if (p == -1)
      * (int *) ptr->value = * (int *) ptr->value ? 0 : 1;
    else
      *(int *) ptr->value = value.SubStr(p + 1).SlowCompare("ON") == 0 ?
	1 : value.SubStr(p + 1).AsInteger();
  else if (ptr->type == LP_DOUBLE_PARAMETER)
    {
      if (p != -1)
	* (double *) ptr->value = value.SubStr(p + 1).AsDouble();
    }
  else if (ptr->type == LP_STRING_PARAMETER)
    {
      if (p != -1)
	* (String *) ptr->value = value.SubStr(p + 1);
    }
    */

bool LongParameters::TranslateExtras(const char * cstr, const char * extras)
{
  if (strchr(cstr, ':') != NULL)
    return false;

  /*
  int option = index.FindStem(cstr);

  if (option == -2)
    {
      // No need to explain ambiguity here ... will be handle by later call
      // to Translate()
      // ExplainAmbiguity(cstr);
      return false;
    }
  */

  StringMapIterator i = index.find(cstr);

  LongParameterList * ptr;

  if ( i != index.end() ) {
    ptr = (LongParameterList *) i->second;
  }
  else {
    i = legacyIndex.find(cstr);
    if ( i == index.end() ) 
      return false;

    ptr = (LongParameterList *) i->second;
    ptr->touched = true;
  }

  if (ptr->type == LP_INT_PARAMETER && CheckInteger(extras)) {
    *(int *) ptr->value = atoi(extras);
    //fprintf(stderr,"*** %s %d\n",extras,*(int*) ptr->value);
    return true;
  }
  else if (ptr->type == LP_DOUBLE_PARAMETER && CheckDouble(extras)) {
    *(double *) ptr->value = atof(extras);
    return true;
  }
  else if (ptr->type == LP_STRING_PARAMETER) {
    *(std::string *) ptr->value = extras;
    return true;
  }

  /*

  if (option >= 0)
    ptr = (LongParameterList *) index.Object(option);
  else
    {
      option = legacyIndex.FindFirstStem(cstr);

      if (option < 0)
	return false;

      ptr = (LongParameterList *) legacyIndex.Object(option);
      ptr->touched = true;
    }

  if (ptr->type == LP_INT_PARAMETER && CheckInteger(extras))
    {
      *(int *) ptr->value = atoi(extras);
      return true;
    }
  else if (ptr->type == LP_DOUBLE_PARAMETER && CheckDouble(extras))
    {
      *(double *) ptr->value = atof(extras);
      return true;
    }
  else if (ptr->type == LP_STRING_PARAMETER)
    {
      *(String *) ptr->value = extras;
      return true;
    }
  */
  return false;
}

void LongParameters::Status(LongParameterList * ptr, int & line_len, bool & need_a_comma)
{
  std::string state;
  int line_start = group_len ? group_len + 5 : 0;

  if (ptr->value == NULL) {
    fprintf(stderr, "%s %*s :", need_a_comma ? "\n" : "", group_len + 2, ptr->description);
    need_a_comma = false;
    line_len = line_start;
  }
  else {
    if (ptr->type == LP_BOOL_PARAMETER)
      state = * (bool *) ptr->value ? " [ON]" : "";
    else if (ptr->type == LP_INT_PARAMETER) {
      //fprintf(stderr,"foo\n");
      if (((* (int *) ptr->value == 1) && (ptr->exclusive)) || (* (int *) ptr->value == 0))
	state = *(int *) ptr->value ? " [ON]" : "";
      else
	catprintf(state," [%d]",*(int*)ptr->value);
      //state = " [", catprintf(state,"* (int *) ptr->value, state += ']';
    }
    else if (ptr->type == LP_DOUBLE_PARAMETER)
      if (* (double *) ptr->value != _NAN_)
	{
	  double value = * (double *) ptr->value;
	  
	  state = " [";
	  if (value == 0.0 || value >= 0.01)
	    catprintf(state, "%.*f", precision, value);
	  else
	    catprintf(state, "%.1e", value);
	  state += ']';
	}
      else
	state = "";
    else if (ptr->type == LP_STRING_PARAMETER)
      state = " [" + * (std::string *) ptr->value + "]";
    
    int item_len = 3 + strlen(ptr->description) + need_a_comma + state.size();
    
    if (item_len + line_len > 78 && line_len > line_start)
	{
	  line_len = line_start;
	  fprintf(stderr, "%s\n%*s", need_a_comma ? "," : "", line_len,  "");
	  need_a_comma = 0;
	  item_len -= 1;
	}

      fprintf(stderr, "%s --%s%s", need_a_comma ? "," : (need_a_comma = true, ""),
	      ptr->description, state.c_str());

      need_a_comma = true;
      line_len += item_len;
    }
}

void LongParameters::Status()
{
  if (description != NULL && description[0] != 0)
    fprintf(stderr, "\n%s\n", description);

  bool need_a_comma = false;
  int  line_len = 0;

  bool legacy_parameters = false;
  int legacy_count = 0;

  for (LongParameterList * ptr = list + 1; ptr->description != NULL; ptr++)
    if (ptr->type == LP_LEGACY_PARAMETERS)
      legacy_parameters = true;
    else if (legacy_parameters == false)
      Status(ptr, line_len, need_a_comma);
    else if (ptr->touched)
      {
	if (legacy_count == 0)
	  {
	    fprintf(stderr, "\n\nAdditional Options:\n %*s ", group_len + 3, "");
	    line_len = group_len + 5;
	    need_a_comma = false;
	  }

	Status(ptr, line_len, need_a_comma);
	legacy_count++;
      }

  fprintf(stderr, "\n");
}

void ParameterList::Add(Parameter * p)
{
  if (count + 1 >= size)
    error("Parameter list size should be increased");

  p->SetWarningBuffer(warnings);
  pl[count++] = p;
};

void ParameterList::Read(int argc, char ** argv, int start)
{
  MakeString(argc, argv, start);
  for (int i=start; i < argc; i++) {
    bool success = false;

    if (argv[i][0] == '-' && argv[i][1])
      for (int j=0; j<count; j++) {
	success = tolower(argv[i][1]) == pl[j]->ch;
	
	if (success) {
	  if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
	    i++;
	  else if (argv[i][2] == 0 && (i+1 < argc) && (argv[i + 1][0] != '-'))
	    pl[j]->Translate(argv[++i]);
	  else
	    pl[j]->Translate(argv[i] + 2);
	  break;
	}
      }
    
    if (!success)
      {
	catprintf(warnings, "Command line parameter %s (#%d) ignored\n", argv[i], i);
      }
  }
}

int ParameterList::ReadWithTrailer(int argc, char ** argv, int start)
{
  MakeString(argc, argv, start);

  int last_success = start - 1;
  bool split = false;

  for (int i=start; i < argc; i++)
    {
      bool success = false;

      if (argv[i][0] == '-' && argv[i][1])
	for (int j=0; j<count; j++)
	  {
	    success = tolower(argv[i][1]) == pl[j]->ch;

	    if (success)
	      {
		if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
		  split = true;
		else if (argv[i][2] == 0 && (i+1 < argc) && (argv[i + 1][0] != '-'))
		  pl[j]->Translate(argv[i + 1]), split = true;
		else
		  pl[j]->Translate(argv[i] + 2);
		break;
	      }
	  }

      if (success)
	for (last_success++; last_success < i; last_success++)
	  catprintf(warnings,"Command line parameter %s (#%d) ignored\n",
			  argv[last_success], last_success);

      if (split)
	{
	  split = false;
	  last_success++;
	  i++;
	}
    }

  return last_success;
};


void ParameterList::Status()
{
  fprintf(stderr, "\nThe following parameters are available.  Ones with \"[]\" are in effect:\n");

  for (int i=0; i<count; i++)
    pl[i]->Status();

  fprintf(stderr, "\n");

  if (warnings.size())
    {
      ::warning("Problems encountered parsing command line:\n\n%s",
		warnings.c_str());
      warnings.clear();
    }

  if (messages.size())
    fprintf(stderr, "NOTES:\n%s\n", messages.c_str());
}

void ParameterList::MakeString(int argc, char ** argv, int start)
{
  int len = 0;

  for (int i=start; i<argc; i++)
    len += strlen(argv[i]) + 1;

  string = new char [len+1];
  string[0] = 0;

  for (int i=start; i<argc; i++)
    {
      strcat(string, argv[i]);
      strcat(string, " ");
    }
}

ParameterList::~ParameterList()
{
  for (int i = 0; i < count; i++)
    delete pl[i];
  delete [] pl;
  delete [] string;
};

bool Parameter::CheckInteger(const char * value)
{
  if (value[0] != '+' && value[0] != '-' &&
      (value[0] < '0' || value[0] > '9'))
    return false;

  int pos = 1;
  while (value[pos] != 0)
    if (value[pos] < '0' || value[pos] > '9')
      return false;
    else
      pos++;

  return true;
}

bool Parameter::CheckDouble(const char * value)
{
  if (value[0] != '+' && value[0] != '-' && value[0] != '.' &&
      (value[0] < '0'  || value[0] > '9'))
    {
      return false;
    }

  bool decimal = value[0] == '.';

  for (int pos = 1; value[pos] != 0; pos++)
    {
      if (value[pos] < '0' || value[pos] > '9')
	{
	  if (!decimal && value[pos] == '.')
	    {
	      decimal = true;
	    }
	  else if (value[pos] == 'e' || value[pos] == 'E')
	    {
	      return CheckInteger(value + pos + 1);
	    }
	}
    }

  return true;
}

void ParameterList::Enforce(bool & var, bool value, const char * format, ...)
{
  if (var == value)
    return;

  var = value;

  va_list ap;
  va_start(ap, format);
  catprintf(messages, format, ap);
  va_end(ap);
}

void ParameterList::Enforce(int & var, int value, const char * format, ...)
{
  if (var == value)
    return;

  var = value;

   va_list ap;
  va_start(ap, format);
  catprintf(messages, format, ap);
  va_end(ap);
}

void ParameterList::Enforce(double & var, double value, const char * format, ...)
{
  if (var == value)
    return;

  var = value;

  va_list ap;
  va_start(ap, format);
  catprintf(messages, format, ap);
  va_end(ap);
}

void ParameterList::Enforce(std::string & var, const char * value, const char * format, ...)
{
  if (var == value)
    return;

  var = value;

  va_list ap;
  va_start(ap, format);
  catprintf(messages, format, ap);
  va_end(ap);
}
